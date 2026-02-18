#' @title Log Rank Test
#' @description Perform a log rank test on a dataframe with support for multiple conditions as separate univariate tests.
#' @param data data.frame, data to perform log rank test on.
#' @param comparisons character, comparisons to perform log rank test on. Can be a single condition or vector of conditions.
#' @param censors character, censors to perform log rank test on.
#' @param censor_prefix character, prefix for censor columns. Default is 'censor_'.
#' @param time_prefix character, prefix for time columns. Default is 'time_to_'.
#' @return data.frame, p-values from log rank test.
#' @importFrom survival survdiff
#' @export
log_rank_test <- function(
    data,
    comparisons,
    censors,
    censor_prefix = "censor_",
    time_prefix = "time_to_") {
  # Validate inputs
  stopifnot(
    "data must be a data.frame" = is.data.frame(data),
    "comparisons must be character vector" = is.character(comparisons),
    "censors must be character vector" = is.character(censors)
  )

  # Check if all comparison variables exist in data
  missing_vars <- comparisons[!comparisons %in% colnames(data)]
  if (length(missing_vars) > 0) {
    stop(sprintf("Comparison variables not found in data: %s", paste(missing_vars, collapse = ", ")))
  }

  # Create an empty dataframe to store the p-values
  p_values <- data.frame(matrix(ncol = length(comparisons), nrow = length(censors)))
  colnames(p_values) <- comparisons
  rownames(p_values) <- censors

  # Loop through each censor
  for (censor in censors) {
    time <- gsub(censor_prefix, time_prefix, censor)

    # Check if censor and time are in the data
    if (!(time %in% colnames(data)) | !(censor %in% colnames(data))) {
      message(sprintf("Missing either time or censor column for %s", censor))
      next
    }

    # Loop through each comparison
    for (i in seq_along(comparisons)) {
      comparison <- comparisons[i]
      # Perform the log rank test
      fmla <- as.formula(sprintf("Surv(%s, %s) ~ %s", time, censor, comparison))
      fit <- do.call(survdiff, list(fmla, data = data))
      p_value <- fit$p

      # Store the p-value
      p_values[censor, comparison] <- p_value
    }
  }

  # Return the p-values dataframe
  return(p_values)
}

#' @title Survival Analysis
#' @description Make survival curves and HR plots for one or more conditions and censors using univariate analyses.
#' @param data data.frame, data to make survival curves from.
#' @param condition character, condition(s) to make survival curves for. Can be a single condition or vector of conditions.
#' @param censors character, censors to make survival curves for.
#' @param outpath character, output directory to save plots to.
#' @param image_type character, type of image to save. Default is 'pdf'.
#' @param censor_prefix character, prefix for censor columns. Default is 'censor_'.
#' @param time_prefix character, prefix for time columns. Default is 'time_to_'.
#' @return list, list of survival plots and HR plots.
#' @importFrom ggsurvfit ggsurvfit add_confidence_interval
#' @importFrom survival Surv survfit survdiff coxph
#' @importFrom ggplot2 coord_cartesian theme_classic theme labs lims element_text ggsave
#' @importFrom broom tidy
#' @importFrom purrr map_dfr
#' @export
survival_analysis <- function(
    data,
    condition,
    censors,
    outpath,
    censor_prefix = "censor_",
    time_prefix = "time_to_") {
  # Validate inputs
  stopifnot(
    "data must be a data.frame" = is.data.frame(data),
    "condition must be character vector" = is.character(condition),
    "censors must be character vector" = is.character(censors)
  )

  # Check if condition variables exist in data
  missing_vars <- condition[!condition %in% colnames(data)]
  if (length(missing_vars) > 0) {
    stop(sprintf("Condition variables not found in data: %s", paste(missing_vars, collapse = ", ")))
  }

  dir.create(file.path(outpath, "survival_plots"), showWarnings = FALSE, recursive = TRUE)
  censors_basenames <- gsub(censor_prefix, "", censors)
  time_vars <- paste0(time_prefix, censors_basenames)

  # Verify the censors and time_vars are in the data
  stopifnot(all(censors %in% colnames(data)))
  stopifnot(all(time_vars %in% colnames(data)))

  # Make lists to store the survival plots and HR results
  survival_plots <- list()
  HR_list <- list()

  # Univariate analysis: analyze each condition separately
  for (cond in condition) {
    for (idx in seq_along(censors)) {
      time_var <- time_vars[idx]
      censor <- censors[idx]
      fmla <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, censor, cond))

      surv_obj <- do.call(survival::survfit, list(fmla, data = data))
      surv_diff <- do.call(survival::survdiff, list(fmla, data = data))
      coxmodel <- do.call(survival::coxph, list(fmla, data = data))

      # Create survival plot
      surv_plot <- ggsurvfit::ggsurvfit(surv_obj) +
        coord_cartesian(clip = "off") +
        ggsurvfit::add_confidence_interval(alpha = 0.1) +
        theme_classic(18) +
        theme(legend.position = "top") +
        lims(x = c(0, max(surv_obj$time) + 100)) +
        labs(
          x = "Time (days)",
          title = sprintf(
            "Survival Analysis: %s ~ %s (n=%d, events=%d)",
            censor, cond, sum(surv_obj$n), sum(surv_obj$n.event)
          )
        ) +
        theme(axis.title.y = element_text(angle = 90, vjust = 0.5))

      # Store plot with appropriate naming
      plot_name <- if (length(condition) == 1) censor else sprintf("%s_%s", cond, censor)
      survival_plots[[plot_name]] <- surv_plot

      # Save plot with appropriate filename
      filename <- if (length(condition) == 1) {
        sprintf("survival_plot_%s.pdf", censor)
      } else {
        sprintf("survival_plot_%s_%s.pdf", cond, censor)
      }
      ggsave(file.path(outpath, "survival_plots", filename), surv_plot)

      # Extract results
      o <- tidy(coxmodel)
      o$censor <- censor
      o$condition_var <- cond
      o$logrank.chisq <- surv_diff$chisq
      o$logrank.pvalue <- surv_diff$pvalue
      o$hazard.ratio <- exp(o$estimate)
      o$ci.upper <- exp(o$estimate + 1.96 * o$std.error)
      o$ci.lower <- exp(o$estimate - 1.96 * o$std.error)
      o <- o[, c("censor", "condition_var", "term", "logrank.chisq", "logrank.pvalue", "hazard.ratio", "ci.lower", "ci.upper", "p.value")]

      HR_list[[plot_name]] <- o
    }
  }

  # Combine HR results
  HR_df <- do.call(rbind, HR_list)
  write.csv(HR_df, file.path(outpath, "hazard_ratios.csv"))

  out <- list(survival_plots = survival_plots, hazard_ratios = HR_df)
  return(out)
}

#' @title Internal Hazard Ratios Helper
#' @description Internal function to calculate hazard ratios for a single condition.
#' @param cs list, list with censor and time variable names.
#' @param data data.frame, data to make survival curves from.
#' @param condition character, a single condition to fit the model for.
#' @param controls character, controls to include in the model. Default is NULL.
#' @return data.frame, hazard ratios for the condition.
#' @importFrom survival coxph cox.zph
#' @importFrom broom tidy
#' @keywords internal
hazards.internal <- function(cs, data, condition, controls) {
  stopifnot("condition must be a single character value" = length(condition) == 1)

  # Clean data - remove NA values right before model fitting
  cols_to_check <- c(condition, cs$censor, cs$time, controls)
  cols_to_check <- cols_to_check[cols_to_check %in% colnames(data)]
  na_count <- sum(!complete.cases(data[, cols_to_check]))
  if (na_count > 0) {
    warning(sprintf("Removing %d rows with NA values", na_count))
    data <- data[complete.cases(data[, cols_to_check]), ]
  }

  fmla <- as.formula(sprintf(
    "survival::Surv(%s, %s) ~ %s",
    cs$time, cs$censor,
    paste0(c(condition, controls), collapse = " + ")
  ))
  coxmodel <- do.call(coxph, list(fmla, data = data))

  # Test proportional hazards assumption
  zph_test <- tryCatch(
    {
      survival::cox.zph(coxmodel)
    },
    error = function(e) {
      warning(sprintf("cox.zph test failed for %s: %s", condition, e$message))
      return(NULL)
    }
  )

  # Extract global test p-value
  zph_pvalue <- if (!is.null(zph_test)) {
    zph_test$table["GLOBAL", "p"]
  } else {
    NA
  }

  # Extract and annotate results
  o <- tidy(coxmodel)
  o <- subset(o, grepl(condition, term))
  o$condition <- condition
  o$censor <- cs$censor
  o$hazard.ratio <- exp(o$estimate)
  o$ci.upper <- exp(o$estimate + 1.96 * o$std.error)
  o$ci.lower <- exp(o$estimate - 1.96 * o$std.error)
  o$n_total <- nrow(data)
  o$n_event <- sum(data[[cs$censor]] == 1)
  o$zph <- zph_pvalue
  o <- o[, c("censor", "condition", "term", "n_total", "n_event", "hazard.ratio", "ci.lower", "ci.upper", "p.value", "zph")]

  return(o)
}

#' @title Hazard Ratios Table
#' @description Make HR tables for one or more conditions and censors using univariate analyses.
#' @param data data.frame, data to make survival curves from.
#' @param condition character, condition(s) to make survival curves for. Can be a single condition or vector of conditions.
#' @param censors character, censors to make survival curves for.
#' @param controls character, controls to include in the model. Default is NULL.
#' @param per_sd logical, whether to standardize the condition. Default is FALSE.
#' @param ovr logical, whether to do one vs rest. Default is FALSE.
#' @param subgroups character vector, variable names for subgroup analyses. Default is NULL.
#' @param censor_prefix character, prefix for censor columns. Default is 'censor_'.
#' @param time_prefix character, prefix for time columns. Default is 'time_to_'.
#' @param verbose logical, whether to print verbose messages. Default is FALSE.
#' @return data.frame, hazard ratios table.
#' @importFrom purrr map_dfr map flatten_dfr
#' @export
hazard_ratios_table <- function(
    data,
    condition,
    censors,
    controls = NULL,
    per_sd = FALSE,
    ovr = FALSE,
    subgroups = NULL,
    censor_prefix = "censor_",
    time_prefix = "time_to_",
    verbose = FALSE) {
  # Validate inputs and prepare data
  stopifnot(
    "data must be a data.frame" = is.data.frame(data),
    "condition must be character vector" = is.character(condition),
    "censors must be character vector" = is.character(censors)
  )

  # Check if condition variables exist in data
  missing_vars <- condition[!condition %in% colnames(data)]
  if (length(missing_vars) > 0) {
    stop(sprintf("Condition variables not found in data: %s", paste(missing_vars, collapse = ", ")))
  }

  censors_basenames <- gsub(censor_prefix, "", censors)
  time_vars <- paste0(time_prefix, censors_basenames)

  stopifnot(
    "Censors not found in data" = all(censors %in% colnames(data)),
    "Time variables not found in data" = all(time_vars %in% colnames(data)),
    "Condition must be numeric for per_sd" = !per_sd || all(sapply(condition, function(x) is.numeric(data[[x]])))
  )

  # Standardize condition if requested
  if (per_sd) {
    for (cond in condition) {
      if (is.numeric(data[[cond]])) {
        data[, cond] <- scale(data[, cond])[, 1]
      }
    }
  }

  # Handle subgroup analysis
  if (!is.null(subgroups)) {
    stopifnot("subgroups must be character vector" = is.character(subgroups))

    if (verbose) message(sprintf("Performing subgroup analysis for: %s", paste(subgroups, collapse = ", ")))

    return(purrr::map_dfr(subgroups, function(subgroup_var) {
      vals <- na.omit(unique(data[[subgroup_var]]))

      purrr::map_dfr(vals, function(val) {
        subset_data <- subset(data, data[[subgroup_var]] == val)

        if (verbose) {
          message(sprintf("  Analyzing %s = %s (N = %d)", subgroup_var, val, nrow(subset_data)))
        }

        tryCatch(
          {
            result <- hazard_ratios_table(
              data = subset_data,
              condition = condition,
              censors = censors,
              controls = controls,
              per_sd = per_sd,
              ovr = ovr,
              subgroups = NULL,
              censor_prefix = censor_prefix,
              time_prefix = time_prefix,
              verbose = verbose
            )
            result$subgroup_var <- subgroup_var
            result$subgroup_val <- val
            return(result)
          },
          error = function(e) {
            if (verbose) message(sprintf("  ERROR for %s = %s: %s", subgroup_var, val, e$message))
            return(data.frame())
          }
        )
      })
    }))
  }

  # Prepare censors data
  censors_df <- data.frame(
    censor = censors,
    time = time_vars,
    stringsAsFactors = FALSE
  )

  # Handle one-vs-rest encoding
  if (ovr) {
    if (verbose) message("Performing one-vs-rest encoding")

    # For multiple conditions with OVR, handle each condition separately
    if (length(condition) > 1) {
      return(purrr::map_dfr(condition, function(cond) {
        data_encoded <- one_hot_encode_ovr(data, cond, binary = FALSE)
        vals <- unique(data[[cond]])
        ovr_columns <- paste0(cond, "_", vals)
        ovr_columns <- ovr_columns[ovr_columns %in% colnames(data_encoded)]

        purrr::map_dfr(ovr_columns, function(ovr_col) {
          hazard_ratios_table(
            data = data_encoded,
            condition = ovr_col,
            censors = censors,
            controls = controls,
            censor_prefix = censor_prefix,
            time_prefix = time_prefix,
            ovr = FALSE,
            verbose = verbose
          )
        })
      }))
    } else {
      # Single condition OVR
      data_encoded <- one_hot_encode_ovr(data, condition[1], binary = FALSE)
      vals <- unique(data[[condition[1]]])
      ovr_columns <- paste0(condition[1], "_", vals)
      ovr_columns <- ovr_columns[ovr_columns %in% colnames(data_encoded)]

      return(purrr::map_dfr(ovr_columns, function(ovr_col) {
        hazard_ratios_table(
          data = data_encoded,
          condition = ovr_col,
          censors = censors,
          controls = controls,
          censor_prefix = censor_prefix,
          time_prefix = time_prefix,
          ovr = FALSE,
          verbose = verbose
        )
      }))
    }
  }

  # Calculate hazard ratios (always univariate)
  if (verbose && length(condition) > 1) message("Performing univariate analysis for multiple conditions")

  # Loop over each censor x condition combination
  purrr::map_dfr(seq_len(nrow(censors_df)), function(i) {
    cs <- list(censor = censors_df$censor[i], time = censors_df$time[i])
    purrr::map_dfr(condition, function(cond) {
      hazards.internal(cs, data, cond, controls)
    })
  })
}

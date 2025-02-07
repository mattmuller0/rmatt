#' @title Survival Analysis Functions
#' @description Functions for survival analysis including log rank tests and hazard ratios.
#' @name survival_functions
#' @importFrom ggsurvfit ggsurvfit add_confidence_interval
#' @importFrom survival Surv survfit survdiff coxph
#' @importFrom ggplot2 ggplot coord_cartesian theme_classic theme labs lims annotate element_text ggsave
#' @importFrom glue glue
#' @importFrom broom tidy
#' @importFrom purrr map
#' @importFrom dplyr filter sym
NULL

#' @title Log Rank Test
#' @description Perform a log rank test on a dataframe.
#' @param data data.frame, data to perform log rank test on.
#' @param comparisons character, comparisons to perform log rank test on.
#' @param censors character, censors to perform log rank test on.
#' @param censor_prefix character, prefix for censor columns. Default is 'C_'.
#' @param time_prefix character, prefix for time columns. Default is 'T_'.
#' @return data.frame, p-values from log rank test.
#' @export
log_rank_test <- function(
    data,
    comparisons,
    censors,
    censor_prefix = "censor_",
    time_prefix = "time_to_") {
  # Create an empty dataframe to store the p-values
  p_values <- data.frame(matrix(ncol = length(comparisons), nrow = length(censors)))
  colnames(p_values) <- comparisons
  rownames(p_values) <- censors

  # Loop through each variable of interest
  for (censor in censors) {
    # Subset the data for the current variable of interest
    time <- gsub(censor_prefix, time_prefix, censor)

    # make sure censor and time are in the data
    if (!(time %in% colnames(data)) | !(censor %in% colnames(data))) {
      message(glue::glue("Missing either time or censor column for {censor}"))
      next
    }

    # Loop through each comparison group
    for (i in seq_along(comparisons)) {
      comparison <- comparisons[i]
      # Perform the log rank test
      fmla <- as.formula(glue::glue("Surv({time}, {censor}) ~ {comparison}"))
      fit <- do.call(survdiff, list(fmla, data = data))
      p_value <- fit$p

      # Store the p-value in the p_values dataframe
      p_values[censor, comparison] <- p_value
    }
  }

  # Return the p-values dataframe
  return(p_values)
}

#' @title Survival Analysis
#' @description Make survival curves and HR plots for a condition and censors.
#' @param df data.frame, data to make survival curves from.
#' @param condition character, condition to make survival curves for.
#' @param censors character, censors to make survival curves for.
#' @param outdir character, output directory to save plots to.
#' @param image_type character, type of image to save. Default is 'pdf'.
#' @param censor_prefix character, prefix for censor columns. Default is 'C_'.
#' @param time_prefix character, prefix for time columns. Default is 'T_'.
#' @return list, list of survival plots and HR plots.
#' @export
survival_analysis <- function(
    df,
    condition,
    censors,
    outdir,
    image_type = "pdf",
    censor_prefix = "censor_",
    time_prefix = "time_to_") {
  dir.create(file.path(outdir, "survival_plots"), showWarnings = FALSE, recursive = TRUE)
  censors_basenames <- gsub(censor_prefix, "", censors)
  time_vars <- paste0(time_prefix, censors_basenames)

  # verify the censors and time_vars are in the data
  stopifnot(all(censors %in% colnames(df)))
  stopifnot(all(time_vars %in% colnames(df)))

  # make a list to store the survival plots
  survival_plots <- list()
  HR_list <- list()

  # loop over the censors
  for (idx in seq_along(censors)) {
    time_var <- time_vars[idx]
    censor <- censors[idx]
    fmla <- as.formula(glue::glue("survival::Surv({time_var}, {censor}) ~ {condition}"))

    surv_obj <- do.call(survival::survfit, list(fmla, data = df))
    surv_diff <- do.call(survival::survdiff, list(fmla, data = df))
    coxmodel <- do.call(survival::coxph, list(fmla, data = df))

    surv_plot <- ggsurvfit::ggsurvfit(surv_obj) +
      coord_cartesian(clip = "off") +
      ggsurvfit::add_confidence_interval(alpha = 0.1) +
      theme_classic(18) +
      theme(legend.position = "top") +
      lims(x = c(0, max(surv_obj$time) + 100)) +
      labs(x = "Time (days)", title = glue::glue("Survival Analysis: {censor} (n={sum(surv_obj$n)})")) +
      theme(axis.title.y = element_text(angle = 90, vjust = 0.5))

    # add the surv_plot to the list
    survival_plots[[censor]] <- surv_plot
    ggplot2::ggsave(filename = file.path(outdir, "survival_plots", glue::glue("survival_plot_{censor}.{image_type}")), plot = surv_plot)

    # make a dataframe of the OR & pvalue and HR & pvalues
    o <- tidy(coxmodel)
    o <- o %>%
      mutate(
        censor = censor,
        logrank.chisq = surv_diff$chisq,
        logrank.pvalue = surv_diff$pvalue,
        hazard.ratio = exp(estimate),
        ci.upper = exp(estimate + 1.96 * std.error),
        ci.lower = exp(estimate - 1.96 * std.error)
      ) %>%
      select(
        censor, term,
        logrank.chisq, logrank.pvalue,
        hazard.ratio, ci.lower, ci.upper, p.value
      )

    HR_list[[censor]] <- o
  }
  HR_df <- do.call(rbind, HR_list)
  write.csv(HR_df, file.path(outdir, "hazard_ratios.csv"))
  out <- list(survival_plots = survival_plots, hazard_ratios = HR_df)
  return(out)
}

#' @title Internal Hazard Ratios Helper
#' @description Internal function to calculate hazard ratios.
#' @param cs list, list of censors and time variables.
hazards.internal <- function(cs) {
  fmla <- as.formula(glue::glue("survival::Surv({cs$time}, {cs$censor}) ~ {paste0(c(condition, controls), collapse = " + ")}"))
  surv_obj <- do.call(survival::survfit, list(fmla, data = df))
  coxmodel <- do.call(survival::coxph, list(fmla, data = df))

  # make a dataframe of the HRs & pvalue and HR & pvalues
  o <- tidy(coxmodel)
  o <- o %>% filter(grepl(condition, term))
  o <- o %>%
    mutate(
      condition = condition,
      censor = cs$censor,
      hazard.ratio = exp(estimate),
      ci.upper = exp(estimate + 1.96 * std.error),
      ci.lower = exp(estimate - 1.96 * std.error)
    ) %>%
    select(
      censor, condition, term,
      hazard.ratio, ci.lower, ci.upper, p.value
    )
  return(o)
}

#' @title Hazard Ratios Table
#' @description Make HR tables for a condition and censors.
#' @param df data.frame, data to make survival curves from.
#' @param condition character, condition to make survival curves for.
#' @param censors character, censors to make survival curves for.
#' @param controls character, controls to include in the model. Default is NULL.
#' @param per_sd logical, whether to standardize the condition. Default is FALSE.
#' @param ovr logical, whether to do one vs rest. Default is FALSE.
#' @param censor_prefix character, prefix for censor columns. Default is 'C_'.
#' @param time_prefix character, prefix for time columns. Default is 'T_'.
#' @return data.frame, hazard ratios table.
#' @importFrom ggsurvfit ggsurvfit
#' @importFrom survival survfit survdiff coxph
#' @importFrom broom tidy
#' @importFrom purrr map
#' @export
hazard_ratios_table <- function(
    df,
    condition,
    censors,
    controls = NULL,
    per_sd = FALSE,
    ovr = FALSE,
    censor_prefix = "censor_",
    time_prefix = "time_to_") {
  censors_basenames <- gsub(censor_prefix, "", censors)
  time_vars <- paste0(time_prefix, censors_basenames)

  # verify the censors and time_vars are in the data
  stopifnot(all(censors %in% colnames(df)))
  stopifnot(all(time_vars %in% colnames(df)))

  # check if there are NA values in the condition
  if (any(is.na(df[, condition]))) {
    warning("NA values in condition, conrols, or censors -- removing")
  }
  df <- df[complete.cases(df[, condition, censors, controls]), ]

  # check if we are per_sd
  if (per_sd) {
    # make sure the condition is continuous
    if (!is.numeric(df[, condition])) {
      stop("condition must be numeric if per_sd is TRUE")
    }
    df[, condition] <- scale(df[, condition])
  }

  # vectorize the operation
  censors_list <- split(data.frame(censor = censors, time = time_vars), seq_along(censors))
  HR_list <- lapply(censors_list, hazards.internal)
  if (!ovr) {
    censors_list <- split(data.frame(censor = censors, time = time_vars), seq_along(censors))
    HR_list <- lapply(censors_list, hazards.internal)
  } else {
    # one hot encode the condition
    message("one hot encoding condition")
    df_ <- one_hot_encode_ovr(df, condition, binary = FALSE)
    vals <- unique(df_[, condition])
    columns_ovr <- colnames(df_)[colnames(df_) %in% paste0(condition, "_", vals)]

    # lapply the function over the columns_ovr and censors
    HR_list <- purrr::map(columns_ovr, ~ hazard_ratios_table(df_, .x, censors, controls = controls, censor_prefix = censor_prefix, time_prefix = time_prefix, ovr = FALSE))
  }

  HR_df <- do.call(rbind, HR_list)
  return(HR_df)
}

#' @title Filtered Hazard Ratio Table
#' @description Make filtered hazard ratio tables.
#' @param data data.frame, data to make survival curves from.
#' @param condition character, condition to make survival curves for.
#' @param risks character, risks to filter by.
#' @param censors character, censors to make survival curves for.
#' @param censor_prefix character, prefix for censor columns. Default is 'C_'.
#' @param time_prefix character, prefix for time columns. Default is 'T_'.
#' @param per_sd logical, whether to standardize the condition. Default is TRUE.
#' @param ovr logical, whether to do one vs rest. Default is FALSE.
#' @param verbose logical, whether to print verbose messages. Default is FALSE.
#' @param ... additional arguments passed to hazard_ratios_table.
#' @return data.frame, filtered hazard ratios table.
#' @importFrom purrr map
#' @importFrom dplyr filter
#' @importFrom glue glue
#' @export
filtered_hazard_ratio_table <- function(
    data,
    condition,
    risks,
    censors,
    censor_prefix = "censor_",
    time_prefix = "time_to_",
    per_sd = TRUE,
    ovr = FALSE,
    verbose = FALSE,
    ...) {
  # make sure the risks are all characters
  if (!all(sapply(risks, is.character))) {
    stop("risks must be characters")
  }

  surv_risk_res <- purrr::map(
    risks,
    function(x) {
      vals <- na.omit(unique(data[, x]))

      if (verbose) {
        message(glue::glue("Filtering for {x}"))
      }

      res <- purrr::map(
        vals,
        function(y) {
          tmp <- dplyr::filter(data, !!sym(x) == y)
          if (verbose) {
            message(glue::glue("    Subfiltering {y}"))
            message(glue::glue("    N = {nrow(tmp)}"))
          }

          tryCatch(
            {
              out <- hazard_ratios_table(
                df = tmp,
                condition = condition,
                censors = censors,
                censor_prefix = censor_prefix,
                time_prefix = time_prefix,
                per_sd = per_sd,
                ovr = ovr,
                ...
              )
              out$x <- x
              out$y <- y
              out$n <- nrow(tmp)
              return(out)
            },
            error = function(e) {
              if (verbose) {
                message(glue::glue("    ERROR: {e}"))
              }
              return(NULL)
            }
          )
        }
      )
      res <- do.call(rbind, res)
    }
  )
  surv_risk_res <- do.call(rbind, surv_risk_res)
  return(surv_risk_res)
}

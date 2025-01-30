#' @title stats_functions
#' @description This script contains functions for statistical analysis.
#' @details The functions include creating a statistics table, hypergeometric scoring, creating a correlation matrix, softmaxing a vector, and min-max normalizing a vector.
#' @name stats_functions
#' @import tableone
#' @import stats
#' @import purrr
NULL

#' Create a Table 1
#' @description Create a stats table from a data.frame
#' @param data data.frame, data to make table 1 from
#' @param groups character, groups to stratify by
#' @param vars character, variables to include in table 1
#' @param printArgs list, arguments to pass to print
#' @param ... other arguments to pass to CreateTableOne
#' @return data.frame, table 1
#' @export
stats_table <- function(
    data,
    groups,
    vars = NULL,
    printArgs = list(
      showAllLevels = FALSE,
      printToggle = FALSE
    ),
    ...) {
  if (!requireNamespace("tableone", quietly = TRUE)) {
    stop("Package 'tableone' is required but not installed.")
  }

  # Get all the variables if none are specified
  if (is.null(vars)) {
    vars <- colnames(data)
    vars <- vars[!vars %in% groups]
  }

  # Make a table 1
  table1 <- tableone::CreateTableOne(
    vars = vars,
    strata = groups,
    data = data,
    addOverall = TRUE,
    ...
  )

  # Make a nice table 1
  table1 <- do.call(print, c(list(table1), printArgs))
  return(table1)
}

#' Hypergeometric Scoring
#' @description Add odds ratio and p-value to a gse object from hypergeometric enrichment analysis
#' @param gse gse object from hypergeometric enrichment analysis in clusterProfiler
#' @param method character, method to use ('fisher' or 'chisq')
#' @param ... other arguments to pass to the test function
#' @return gse object with odds ratio and p-value added
hypergeometric_scoring <- function(
    gse,
    method = "fisher",
    ...) {
  # Get the geneset
  geneset <- gse@gene

  # Get the enrichment sets
  enrichment <- gse@geneSets

  # Get the universe genes
  universe <- gse@universe

  # Make a dataframe where each row is the universe
  df <- data.frame(
    row.names = universe,
    in_geneset = universe %in% geneset
  )
  # Now add each enrichment set as a new column
  for (i in 1:length(enrichment)) {
    df[, names(enrichment)[i]] <- universe %in% enrichment[[i]]
  }

  # Error handling on the method
  if (!method %in% c("fisher", "chisq")) {
    stop("method must be either fisher or chisq")
  }

  # Let's do the OR test for each enrichment set
  if (method == "fisher") {
    ORs <- sapply(
      df[, -1],
      function(x) stats::fisher.test(table(x, df[, "in_geneset"]), ...)$estimate
    )
  } else if (method == "chisq") {
    ORs <- sapply(
      df[, -1],
      function(x) stats::chisq.test(table(x, df[, "in_geneset"]), ...)$estimate
    )
  }

  # Fix the names
  names(ORs) <- gsub("\\..*", "", names(ORs))

  # Add the ORs to the enrichment results
  results_pathways <- rownames(gse@result)
  gse@result$odds_ratio <- ORs[results_pathways]

  # Return the gse object
  return(gse)
}

#' Create a Correlation Matrix
#' @description Create a correlation matrix from a data.frame
#' @param data data.frame, data to make correlation matrix from
#' @param vars1 character, variables to correlate
#' @param vars2 character, variables to correlate
#' @param method character, correlation method to use
#' @param use character, method for handling missing values
#' @param ... other arguments to pass to cor.test
#' @return list with correlation matrix, p-value matrix, and filtered correlation matrix
#' @export
correlation_matrix <- function(data, vars1, vars2, method = "pearson", use = "pairwise.complete.obs", ...) {
  # Make a placeholder for the correlation matrix
  cor_mat <- data.frame()
  p_mat <- data.frame()
  # Map the correlation function over the data
  cor_res <- purrr::map(
    vars1,
    function(x) {
      purrr::map(
        vars2,
        function(y) {
          res <- stats::cor.test(data[, x], data[, y], method = method, use = use, ...)
          cor_mat[x, y] <<- res$estimate
          p_mat[x, y] <<- res$p.value
        }
      )
    }
  )
  # Make a correlation matrix filtered by p-value
  cor_mat_filtr <- cor_mat
  cor_mat_filtr[p_mat > 0.05] <- NA
  return(list(cor_matr = cor_mat, pvalue_matr = p_mat, cor_matr_filtr = cor_mat_filtr))
}

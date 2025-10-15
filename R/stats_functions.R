#' Create a Table 1
#' @description Create a comprehensive descriptive statistics table from a data.frame
#' @param data data.frame, input data for statistical analysis
#' @param groups character, column name(s) for group stratification
#' @param vars character vector, optional list of variables to include in analysis (if NULL, uses all columns except groups)
#' @param table.includeNA logical, whether to include NA as a valid category
#' @param table.addOverall logical, whether to add an overall column
#' @param table.test logical, whether to include statistical tests
#' @param print.showAllLevels logical, whether to show all levels of categorical variables
#' @param print.printToggle logical, whether to show variable labels
#' @param print.catDigits integer, number of digits for categorical variables
#' @param print.contDigits integer, number of digits for continuous variables
#' @param print.pDigits integer, number of digits for p-values
#' @param print.missing logical, whether to show missing data information
#' @param print.nonnormal character vector, variables to be summarized with median/IQR
#' @return matrix containing descriptive statistics table
#' 
#' @export
stats_table <- function(
  data,
  groups,
  vars = NULL,
  # Some parameters to pass to CreateTableOne
  table.includeNA = FALSE,
  table.addOverall = TRUE,
  table.test = TRUE,
  # Some arguments to pass to print
  print.showAllLevels = TRUE,
  print.printToggle = FALSE,
  print.catDigits = 1,
  print.contDigits = 2,
  print.pDigits = 3,
  print.missing = FALSE,
  print.nonnormal = NULL) {
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
  addOverall = table.addOverall,
  test = table.test,
  includeNA = table.includeNA
  )

  # Make a nice table 1
  table1 <- print(
  table1,
  showAllLevels = print.showAllLevels,
  printToggle = print.printToggle,
  catDigits = print.catDigits,
  contDigits = print.contDigits,
  pDigits = print.pDigits,
  missing = print.missing,
  nonnormal = print.nonnormal
  )
  return(table1)
}

#' Hypergeometric Scoring
#' @description Add odds ratio and p-value to a gse object from hypergeometric enrichment analysis
#' @param gse gse object from hypergeometric enrichment analysis in clusterProfiler
#' @param method character, method to use ('fisher' or 'chisq')
#' @param ... other arguments to pass to the test function
#' @return gse object with odds ratio and p-value added
#' @export
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
#' @importFrom purrr map
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
      res <- cor.test(data[, x], data[, y], method = method, use = use, ...)
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

#' Van der Waerden Inverse Normal Transformation
#' @description Applies the van der Waerden inverse normal transformation to a numeric vector.
#' @param x numeric vector to transform
#' @param c numeric constant to add to the denominator (default is 0.375)
#' @return numeric vector with transformed values
#' @export
inverse_normal_transform <- function(x, c = 0) {
  if (!is.numeric(x)) stop("Input x must be numeric.")
  n <- sum(!is.na(x))
  ranks <- rank(x, na.last = "keep", ties.method = "average")
  transformed <- qnorm((ranks - c) / (n + 1 - 2 * c))
  return(transformed)
}
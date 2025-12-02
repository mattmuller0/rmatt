#' Function to get genes from a results table
#' @param results results table
#' @param pval p-value cutoff
#' @param log2fc log2 fold change cutoff
#' @param gene_col gene column name
#' @param pval_col p-value column name
#' @param log2fc_col log2 fold change column name
#' @return dataframe of genes separated by up and down
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter pull
#' @export
getGenes <- function(results, pval = 0.05, metric = 0, name_col = "rownames", pval_col = "padj", metric_col = "log2FoldChange") {
  if (name_col == "rownames") {
  results <- rownames_to_column(results, var = name_col)
  }
  up <- results %>%
  filter(!!sym(pval_col) < pval & !!sym(metric_col) > metric) %>%
  pull(!!sym(name_col))
  down <- results %>%
  filter(!!sym(pval_col) < pval & !!sym(metric_col) < -metric) %>%
  pull(!!sym(name_col))
  out <- data.frame(features = c(up, down), direction = c(rep("up", length(up)), rep("down", length(down))))
  return(out)
}

#' Function to add missing rows to a matrix
#' @param data matrix, rows = genes, cols = samples
#' @param rows vector of rownames to add
#' @param sorted logical, whether to sort the rows alphabetically
#' @return data matrix, rows = genes, cols = samples
#' @export
add_missing_rows <- function(
  data,
  rows,
  sorted = TRUE) {
  missingRowNames <- rows[which(!rows %in% rownames(data))]
  data_tmp <- as.data.frame(matrix(0, nrow = length(missingRowNames), ncol = dim(data)[2]))
  colnames(data_tmp) <- colnames(data)
  rownames(data_tmp) <- missingRowNames
  data <- rbind(data_tmp, data)
  if (sorted) {
  data <- data[order(rownames(data)), ]
  }
  return(data)
}

#' Function to make a SummarizedExperiment object
#' @param countsMatr matrix of counts, rows = genes, cols = samples
#' @param colData data.frame of sample metadata, rows = samples
#' @return SummarizedExperiment object
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
make_se <- function(countsMatr, colData) {
  sample_ids <- intersect(colnames(countsMatr), row.names(colData))
  se <- SummarizedExperiment(assays = list(counts = as.matrix(countsMatr[, sample_ids])), colData = colData[sample_ids, ])
  return(se)
}

#' Function to make a DESeqDataSet object
#' @param countsMatr matrix of counts, rows = genes, cols = samples
#' @param colData data.frame of sample metadata, rows = samples
#' @param design formula for the design matrix
#' @return DESeqDataSet object
#' 
#' @export
make_dds <- function(countsMatr, colData, design = ~ 1) {
  sample_ids <- intersect(colnames(countsMatr), row.names(colData))
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(countsMatr[, sample_ids]), colData = colData[sample_ids, ], design = design)
  return(dds)
}

#' Function to save a SummarizedExperiment object into multiple files if the slots are filled in the object
#' @param se SummarizedExperiment object
#' @param path path to save files
#' @param normalize how to normalize the counts
#' @importFrom SummarizedExperiment colData rowData
#' @export
save_se <- function(se, path, normalize = "mor") {
  dir.create(path, showWarnings = F, recursive = T)

  # Normalize counts based on the specified method
  counts <- normalize_counts(se, method = normalize)
  colData <- as.data.frame(colData(se))
  rowData <- as.data.frame(rowData(se))

  # Save the counts, colData, and rowData to CSV files
  write.csv(counts, file = file.path(path, paste0("counts_", normalize, ".csv")), quote = F, row.names = T)
  if (!is.null(colData)) {
  write.csv(colData, file = file.path(path, "metadata.csv"), quote = T, row.names = T)
  }
  if (!is.null(rowData)) {
  write.csv(rowData, file = file.path(path, "rowData.csv"), quote = T, row.names = T)
  }
}

#' Function to remove NA variables from a summarized experiment object
#' @param se SummarizedExperiment object
#' @param columns list of columns to remove NAs from
#' @return DESeqDataSet object with NAs removed
#' @importFrom SummarizedExperiment colData assay
#' @importFrom S4Vectors DataFrame
#' @importFrom tidyr drop_na any_of
remove_na_variables <- function(se, columns) {
  col_data <- as.data.frame(colData(se))
  assay_data <- assay(se)
  col_data <- drop_na(col_data, any_of(columns))
  col_data <- DataFrame(col_data)
  new_se <- make_se(assay_data, col_data)
  return(new_se)
}

#' Function to get the number of samples in each group
#' @param metadata dataframe, metadata
#' @param id character, column name of the sample ID
#' @param events_term character, prefix of the events columns
#' @param subset character, column name of the grouping variable
#' @return dataframe with the number of samples in each group
#' @importFrom dplyr select mutate summarise_all starts_with
get_events_breakdown <- function(metadata, id = "PATNUM", events_term = "C_", subset = NULL) {
  metadata_patnums <- metadata[, id]
  if (!is.null(subset)) {
  patnum_match <- metadata_patnums[metadata[, subset] == 1]
  } else {
  patnum_match <- metadata_patnums
  }
  if (is.na(patnum_match[1])) {
  stop("No matches found from the PATNUMs provided")
  }
  subtype_breakdown <- metadata[patnum_match, ] %>%
  select(starts_with(events_term)) %>%
  summarise_all(~ sum(., na.rm = TRUE)) %>%
  mutate(total = rowSums(.))
  return(subtype_breakdown)
}

#' Function to get make pairwise combinations
#' @param vec vector of values
#' @return list of pairwise combinations
#' @export
pairwise_combos <- function(vec) {
  unique_vals <- unique(vec)
  combos <- combn(unique_vals, 2)
  list_of_combos <- list()
  for (i in 1:ncol(combos)) {
  list_of_combos[[i]] <- combos[, i]
  }
  return(list_of_combos)
}

#' Function to one hot encode dataframe column
#' @param data dataframe
#' @param column character, column name
#' @param binary logical, whether to make the column binary
#' @return dataframe with one hot encoded column
#' @importFrom glue glue
#' @export
one_hot_encode_ovr <- function(data, column, binary = TRUE) {
  unique_vals <- unique(data[[column]])
  for (val in unique_vals) {
  new_col_name <- paste0(column, "_", val)
  if (binary) {
    data[[new_col_name]] <- ifelse(data[[column]] == val, 1, 0)
  } else {
    data[[new_col_name]] <- factor(ifelse(data[[column]] == val, val, "rest"), levels = c("rest", val))
  }
  }
  return(data)
}

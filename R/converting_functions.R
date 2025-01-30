#' @title Converting Functions
#' @description This script contains functions for converting gene IDs and normalizing data.
#' @details The functions include mapping gene IDs using `AnnotationDbi`, converting gene IDs using `biomaRt`, detecting gene ID types, and normalizing data.
#' @name converting_functions
#' @importFrom org.Hs.eg.db org.Hs.eg.db 
#' @importFrom AnnotationDbi mapIds
#' @importFrom DESeq2 sizeFactors estimateSizeFactors counts varianceStabilizingTransformation rlog
#' @importFrom SummarizedExperiment assay
#' @importFrom edgeR cpm rpkm calcNormFactors
#' @importFrom singscore rankGenes
#' @importFrom stringr str_split
#' @importFrom glue glue
#' @importFrom dplyr %>%
NULL

#' Map Gene IDs
#'
#' @description Function to map gene IDs with `AnnotationDbi`.
#' @param geneList list or vector of gene IDs.
#' @param from character, type of gene ID to convert from (if NULL, will detect).
#' @param to character, type of gene ID to convert to.
#' @param orgDb OrgDb object, organism database to use (default is `org.Hs.eg.db`).
#' @param remove_missing logical, keep unmatched gene IDs.
#' @param ... additional arguments to pass to `AnnotationDbi::mapIds`.
#' @return list, converted gene IDs.
map_gene_ids <- function(geneList, from = NULL, to, orgDb = org.Hs.eg.db, remove_missing = FALSE, ...) {
  if (is.null(from)) {
    from <- detect_gene_id_type(geneList)
  }

  geneList_copy <- geneList

  if (grepl("\\.", geneList[1])) {
    geneList <- sapply(geneList, function(x) stringr::str_split(x, "\\.", simplify = TRUE)[1])
  }

  out_list <- AnnotationDbi::mapIds(orgDb, keys = geneList, column = to, keytype = from, ...)

  if (remove_missing) {
    out_list <- out_list[!is.na(out_list)]
  } else {
    out_list[is.na(out_list)] <- geneList_copy[is.na(out_list)] %>% as.vector()
  }

  return(out_list)
}

#' Detect Gene ID Type
#'
#' @description Function to detect gene ID type based on the first gene ID.
#' @param geneList list or vector of gene IDs.
#' @param strip logical, whether to strip version numbers from gene IDs.
#' @return character, type of gene ID.
detect_gene_id_type <- function(geneList, strip = TRUE) {
  id_types <- list(
    ENSEMBL = "^ENSG[0-9]+$",
    ENSEMBLTRANS = "^ENST[0-9]+$",
    ENSEMBLPROT = "^ENSP[0-9]+$",
    ENTREZID = "^[0-9]+$",
    IMAGE = "^IMAGE:[0-9]+$",
    GOID = "^GO:[0-9]+$",
    PFAM = "^PF[0-9]+$",
    REFSEQ = "^N[MP]_[0-9]+$",
    ENZYME = "^[0-9]+(\\.(([0-9]+)|-)+)3$",
    MAP = "^[0-9XY]+((([pq])|(cen))(([0-9]+(\\.[0-9]+)?)|(ter))?(-([0-9XY]+)?(([pq]?)|(cen))((ter)|([0-9]+(\\.[0-9]+)?))?)?)?$",
    UNIGENE = "^Hs\\.([0-9]+)$",
    SYMBOL = "^[A-Z][A-Z0-9]*(_[A-Z0-9]+)*(-[A-Z0-9]+)*$",
    GENEBANK_Nucleotide = "^[A-Z][0-9]5$",
    GENEBANK_Protein = "^[A-Z]3[0-9]5$",
    GENEBANK_WGS = "^[A-Z]4[0-9]8[0-9]?[0-9]?$",
    GENEBANK_MGA = "^[A-Z]5[0-9]7$",
    GENENAME = " ",
    .Affymetrix = "(^AFFX-)|(^[0-9]+_([abfgilrsx]_)?([as]t)|(i))$",
    .Illumina = "^ILMN_[0-9]+$",
    .Agilent = "^A_[0-9]+_P[0-9]+$"
  )

  if (grepl("\\.", geneList[1]) && strip) {
    geneList <- sapply(geneList, function(x) stringr::str_split(x, "\\.", simplify = TRUE)[1])
    message("Gene ID version numbers have been stripped.")
  }

  matches <- sapply(id_types, function(x) grepl(x, geneList))

  type <- colnames(matches)[which(matches[1, ] == TRUE)]

  if (length(type) > 1) {
    warning(glue::glue("Multiple gene ID types detected. Returning first match: {type[1]}"))
    type <- type[1]
  }

  if (length(type) == 0) {
    stop("No gene ID type detected.")
  }

  return(type)
}

#' Function to normalize dds object and return counts
#'
#' @param dds DESeq2 object
#' @param method normalization method
#' @param log2 logical, whether to log2 transform the counts
#' @return normalized counts
#' @export
normalize_counts <- function(dds, method = "mor", log2 = FALSE) {
  options <- c("mor", "vst", "vsd", "log2", "rld", "cpm", "rlog", "rpkm", "none", "tmm", "log2-mor", "rank")
  if (method %in% options) {
    normalize <- method
  } else {
    stop("Invalid normalization method. Please choose from: log2-mor, mor, vst, vsd, log2, rld, cpm, rlog, rpkm, tmm, rank, none")
  }
  if (normalize == "mor") {
    if (is.null(DESeq2::sizeFactors(dds))) {
      dds <- DESeq2::estimateSizeFactors(dds)
    }
    counttable <- DESeq2::counts(dds, normalize = TRUE)
  }
  if (normalize %in% c("vsd", "vst")) {
    counttable <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds))
  }
  if (normalize == "log2") {
    counttable <- log2(SummarizedExperiment::assay(dds) + 1)
  }
  if (normalize == "log2-mor") {
    counttable <- log2(DESeq2::counts(dds, normalize = TRUE) + 1)
  }
  if (normalize %in% c("rld", "rlog")) {
    counttable <- DESeq2::rlog(dds)
  }
  if (normalize == "cpm") {
    counttable <- edgeR::cpm(dds)
  }
  if (normalize == "rpkm") {
    counttable <- edgeR::rpkm(dds)
  }
  if (normalize == "tmm") {
    counttable <- edgeR::cpm(edgeR::calcNormFactors(dds, method = "TMM"))
  }
  if (normalize == "rank") {
    counttable <- singscore::rankGenes(dds)
  }
  if (normalize == "none") {
    counttable <- SummarizedExperiment::assay(dds)
  }
  counts <- as.data.frame(counttable)
  if (log2) {
    counts <- log2(counts + 1)
  }
  return(counts)
}

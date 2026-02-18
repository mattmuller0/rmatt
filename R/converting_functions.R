#' Map Gene IDs
#' @description Function to map gene IDs with `AnnotationDbi`.
#' @param geneList list or vector of gene IDs.
#' @param from character, type of gene ID to convert from (if NULL, will detect).
#' @param to character, type of gene ID to convert to.
#' @param orgDb OrgDb object, organism database to use (default is `org.Hs.eg.db`).
#' @param remove_missing logical, keep unmatched gene IDs.
#' @param ... additional arguments to pass to `AnnotationDbi::mapIds`.
#' @return list, converted gene IDs.
#' @importFrom stringr str_split
#' @export
map_gene_ids <- function(geneList, from = NULL, to, orgDb = NULL, remove_missing = FALSE, ...) {
  # Check for org.Hs.eg.db if no orgDb provided
  if (is.null(orgDb)) {
    check_bioc_package("org.Hs.eg.db")
    orgDb <- org.Hs.eg.db::org.Hs.eg.db
  }
  
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
    out_list[is.na(out_list)] <- as.vector(geneList_copy[is.na(out_list)])
  }

  return(out_list)
}

#' Detect Gene ID Type
#' @description Function to detect gene ID type based on the first gene ID.
#' @param geneList list or vector of gene IDs.
#' @param strip logical, whether to strip version numbers from gene IDs.
#' @return character, type of gene ID.
#' @importFrom stringr str_split
#' @importFrom glue glue
#' @export
detect_gene_id_type <- function(geneList, strip = TRUE) {
  id_types <- list(
    # Start with the most common gene ID types
    ENSEMBL = "^ENSG[0-9]+$",
    ENSEMBLTRANS = "^ENST[0-9]+$",
    ENSEMBLPROT = "^ENSP[0-9]+$",
    SYMBOL = "^[A-Z][A-Z0-9]*(_[A-Z0-9]+)*(-[A-Z0-9]+)*$",
    ENTREZID = "^[0-9]+$",
    REFSEQ = "^N[MP]_[0-9]+$",
    UNIGENE = "^Hs\\.([0-9]+)$",
    # Now some less common gene ID types
    IMAGE = "^IMAGE:[0-9]+$",
    GOID = "^GO:[0-9]+$",
    PFAM = "^PF[0-9]+$",
    ENZYME = "^[0-9]+(\\.(([0-9]+)|-)+)3$",
    MAP = "^[0-9XY]+((([pq])|(cen))(([0-9]+(\\.[0-9]+)?)|(ter))?(-([0-9XY]+)?(([pq]?)|(cen))((ter)|([0-9]+(\\.[0-9]+)?))?)?)?$",
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
    warning("Gene ID version numbers have been stripped.")
  }

  matches <- sapply(id_types, function(x) grepl(x, geneList, perl = TRUE))
  
  # get the colSums of the matches
  matched_num <- colSums(matches)
  type <- names(id_types)[which.max(matched_num)]

  if (length(type) > 1) {
    message(glue::glue("Multiple gene ID types detected. Returning first match: {type[1]} from ({paste(type, collapse = ', ')})"))
    type <- type[1]
  }

  if (all(matched_num == 0)) {
    stop("No gene ID type detected.")
  }

  return(type)
}

#' Fetch Gene Coordinates from Bioconductor
#' @description Function to fetch gene coordinates from Bioconductor based on gene symbols and genome build.
#' @param genes Character vector of gene symbols
#' @param genome_build Genome build (default: "hg38")
#' @return Data frame with gene, chr, start, end columns
#' @importFrom dplyr filter distinct select inner_join bind_rows
#' @importFrom AnnotationDbi select
#' @importFrom GenomicFeatures genes
#' @importFrom glue glue
#' @export
fetch_gene_coordinates <- function(genes, key = "SYMBOL", genome_build = "hg38", filter_chromosomes = "^chr([1-9]|1[0-9]|2[0-2]|X|Y)$") {
    # Check for required packages
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop("Package 'org.Hs.eg.db' is required. Install with BiocManager::install('org.Hs.eg.db')")
    }
    
    txdb_pkg <- switch(
        genome_build,
        "hg38" = "TxDb.Hsapiens.UCSC.hg38.knownGene",
        "hg19" = "TxDb.Hsapiens.UCSC.hg19.knownGene",
        stop(glue("Unsupported genome build: {genome_build}"))
    )
    
    if (!requireNamespace(txdb_pkg, quietly = TRUE)) {
        stop(glue("Package '{txdb_pkg}' is required. Install with BiocManager::install('{txdb_pkg}')"))
    }
    
    message(glue("Fetching coordinates for {length(genes)} genes using {genome_build}..."))
    
    # Map symbols to Entrez IDs
    key_to_entrez <- AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = genes,
        columns = c("ENTREZID", key),
        keytype = key
    ) %>%
        filter(!is.na(ENTREZID)) %>%
        distinct(ENTREZID, .keep_all = TRUE)
    
    # Get gene ranges from TxDb
    txdb <- get(txdb_pkg, envir = asNamespace(txdb_pkg))
    gene_ranges <- GenomicFeatures::genes(txdb)
    
    # Filter to our genes
    matched_ranges <- gene_ranges[names(gene_ranges) %in% symbol_to_entrez$ENTREZID]
    
    # Convert to data frame
    coords_df <- as.data.frame(matched_ranges) %>%
        rownames_to_column("ENTREZID") %>%
        select(ENTREZID, chr = seqnames, start, end) %>%
        inner_join(key_to_entrez, by = "ENTREZID") %>%
        select(gene = !!sym(key), chr, start, end)

    # Filter to standard chromosomes
    chromosomes <- unique(as.character(coords_df$chr))
    standard_chromosomes <- chromosomes[grepl(filter_chromosomes, chromosomes)]
    if (length(standard_chromosomes) == 0) {
        warning(glue("No standard chromosomes found for genome build '{genome_build}'"))
        standard_chroms <- chromosomes
    }

    # Filter to standard chromosomes
    coords_df <- coords_df %>%
        filter(as.character(chr) %in% standard_chromosomes)
    
    # Add missing genes as NA
    missing <- setdiff(genes, coords_df$gene)
    if (length(missing) > 0) {
        missing_df <- data.frame(
            gene = missing,
            chr = NA_character_,
            start = NA_integer_,
            end = NA_integer_
        )
        coords_df <- bind_rows(coords_df, missing_df)
    }

    # Make the chr a factor with standard order
    coords_df$chr <- factor(coords_df$chr, levels = standard_chromosomes)
    
    return(coords_df)
}

#' Function to normalize dds object and return counts
#' @description Function to normalize DESeq2 dds object and return counts.
#' @param dds DESeq2 or SummarizedExperiment object
#' @param method Normalization method. One of: "mor" (median of ratios), "vst"/"vsd" 
#'   (variance stabilizing transformation), "log2", "log2-mor", "rld"/"rlog" 
#'   (regularized log), "cpm" (counts per million), "rpkm", "tmm", "rank", or "none"
#' @param log2 Logical, whether to additionally log2 transform the counts (default: FALSE)
#' @return Data frame of normalized counts with genes as rows and samples as columns
#' @importFrom SummarizedExperiment assay
#' @export
normalize_counts <- function(
  dds, 
  method = c("mor", "vst", "vsd", "log2", "log2-mor", "rld", "rlog", "cpm", "rpkm", "tmm", "rank", "none"), 
  log2 = FALSE) {
  method <- match.arg(method)
  
  # Handle vsd as alias for vst, rld as alias for rlog
  method <- switch(method,
    vsd = "vst",
    rld = "rlog",
    method
  )
  
  counttable <- switch(method,
    mor = {
      if (is.null(DESeq2::sizeFactors(dds))) {
        dds <- DESeq2::estimateSizeFactors(dds)
      }
      DESeq2::counts(dds, normalize = TRUE)
    },
    vst = SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds)),
    log2 = log2(SummarizedExperiment::assay(dds) + 1),
    `log2-mor` = {
      if (is.null(DESeq2::sizeFactors(dds))) {
        dds <- DESeq2::estimateSizeFactors(dds)
      }
      log2(DESeq2::counts(dds, normalize = TRUE) + 1)
    },
    rlog = SummarizedExperiment::assay(DESeq2::rlog(dds)),
    cpm = edgeR::cpm(dds),
    rpkm = edgeR::rpkm(dds),
    tmm = edgeR::cpm(edgeR::calcNormFactors(dds, method = "TMM")),
    rank = singscore::rankGenes(dds),
    none = SummarizedExperiment::assay(dds)
  )
  
  counts <- as.data.frame(counttable)
  if (log2) {
    counts <- log2(counts + 1)
  }
  return(counts)
}


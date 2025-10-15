# Some general exports I tend to need
# Note: Re-exporting these functions forces SummarizedExperiment to load on package startup.
# Users can call SummarizedExperiment::SummarizedExperiment() directly instead.
# These exports have been removed to improve package load time.

# #' @export
# SummarizedExperiment::SummarizedExperiment
#
# #' @export
# SummarizedExperiment::assay
#
# #' @export
# SummarizedExperiment::rowData
#
# #' @export
# SummarizedExperiment::colData

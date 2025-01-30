#' @title Phenomics Functions
#' @description Functions for EHR / Biobank data processing.
#' @name phenomics_functions
#' @importFrom dplyr filter pull mutate group_by summarise
#' @importFrom glue glue
#' @importFrom purrr map
#' @importFrom tidyr pivot_longer pivot_wider
NULL

#' @title Make Composite Coding CSV
#' @description Function to take any codes and make a composite variable CSV.
#' @param lol List of lists containing the codes.
#' @param composite_name Name of the composite variable.
#' @return A data frame with the composite variable.
#' @export
make_compsite_coding_csv <- function(lol, composite_name) {
    # ensure the lol args are all named vectors
    if (!all(sapply(list(lol), is.vector))) {
        stop("All unnamed arguments must be named vectors")
    }
    out <- data.frame(
        key = rep(composite_name, length(unlist(unlist(lol)))),
        field = c(rep(names(lol), sapply(lol, length))),
        value = c(unlist(unlist(lol)))
    )
    return(out)
}

#' @title Make Phenotypes
#' @description Function to make a phenotype from the composite encoding.
#' @param hesin Data frame of hesin data pivoted longer [dnx_hesin_id, Participant ID, field, value].
#' @param encoding Data frame of composite encoding [key, field, value].
#' @return A data frame with composite encoding.
#' @export
make_phenotypes <- function(hesin, encoding) {
    # ensure the required columns are present
    reqs <- c("Participant ID", "field", "value")
    if (!all(reqs %in% colnames(hesin))) {
        stop(glue::glue("Missing required columns in hesin data\n[{paste0(reqs, collapse = ', ')}]"))
    }
    reqs <- c("key", "field", "value")
    if (!all(reqs %in% colnames(encoding))) {
        stop(glue::glue("Missing required columns in encoding data\n[{paste0(reqs, collapse = ', ')}]"))
    }

    out <- list()
    for (i in seq_along(unique(encoding$key))) {
        k <- unique(encoding$key)[i]
        print(glue("Processing {unique(encoding$key)[i]}"))
        tmp <- encoding %>% filter(key == k)
        fields <- unique(tmp$field)
        o <- list()
        for (j in seq_along(fields)) {
            print(glue("    Field: {fields[j]}"))
            v <- unique(tmp %>% filter(field == fields[j]) %>% pull(value))
            v_reg <- paste0("^(", paste(v, collapse = "|"), ")")
            print(glue("        Values: {paste0(v, collapse = '|')}"))
            o[[j]] <- hesin %>%
                filter(field == fields[j]) %>%
                filter(grepl(v_reg, value)) %>%
                mutate(phenotype = k)
            print(glue("        N: {nrow(o[[j]])}"))
        }
        out[[i]] <- do.call(rbind, o)
    }
    out <- do.call(rbind, out)
    summary <- out %>%
        group_by(phenotype) %>%
        summarise(n = n())
    print("Summary of phenotypes")
    print(summary)
    return(out)
}

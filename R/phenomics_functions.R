#' @title Phenomics Functions
#' @description Functions for EHR / Biobank data processing.
#' @name phenomics_functions
#' @importFrom dplyr filter pull mutate group_by summarise n_distinct bind_rows
#' @importFrom glue glue
#' @importFrom purrr map
#' @importFrom tidyr pivot_longer pivot_wider
NULL

#' Make Composite Coding CSV
#' @param lol List of lists containing the codes
#' @param composite_name Name of the composite variable
#' @return Data frame with composite variable mapping
#' @export
make_composite_coding <- function(lol, composite_name) {
  if (!is.list(lol) || !all(sapply(lol, is.vector))) {
    stop("Input must be a list of named vectors")
  }
  
  if (length(lol) == 0) {
    stop("Input list is empty")
  }
  
  results <- map(names(lol), function(field) {
    data.frame(
      key = composite_name,
      field = field,
      value = unlist(lol[[field]]),
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(results)
}

#' Make Phenotypes from Composite Encoding
#' @param hesin Data frame with columns [dnx_hesin_id, Participant ID, field, value]
#' @param encoding Data frame with columns [key, field, value]
#' @return Data frame with matched phenotypes
#' @export
make_phenotypes <- function(hesin, encoding) {
  required_hesin <- c("field", "value")
  required_encoding <- c("key", "field", "value")
  
  if (!all(required_hesin %in% colnames(hesin))) {
    stop("Missing hesin columns: ", 
         paste(setdiff(required_hesin, colnames(hesin)), collapse = ", "))
  }
  
  if (!all(required_encoding %in% colnames(encoding))) {
    stop("Missing encoding columns: ", 
         paste(setdiff(required_encoding, colnames(encoding)), collapse = ", "))
  }
  
  phenotype_results <- map(unique(encoding$key), function(phenotype) {
    message(glue("Processing {phenotype}"))
    phenotype_encoding <- encoding[encoding$key == phenotype, ]
    
    field_results <- map(unique(phenotype_encoding$field), function(field) {
      values <- unique(phenotype_encoding$value[phenotype_encoding$field == field])
      value_pattern <- paste0("^(", paste(values, collapse = "|"), ")")
      
      matches <- hesin %>%
        filter(field == !!field, grepl(value_pattern, value)) %>%
        mutate(phenotype = phenotype)
      
      message(glue("  {field}: {n_distinct(matches$`Participant ID`)} participants"))
      matches
    })
    
    bind_rows(field_results)
  })
  
  bind_rows(phenotype_results) %>%
    group_by(phenotype) %>%
    mutate(n_codes = n_distinct(value)) %>%
    ungroup()
}
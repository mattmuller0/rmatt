#' Install all required packages from a script
#' @description Install all required packages from a script
#' @param script character, path to script
#' @return none
#' @importFrom BiocManager install
#' @export
install_packages_from_script <- function(script) {
  script <- utils::readLines(script)

  packages <- script[grepl("library", script)]
  packages <- gsub("library\\(", "", packages)
  packages <- gsub("suppressMessages\\(", "", packages)
  packages <- gsub("\\)", "", packages)
  packages <- gsub('"', "", packages)
  packages <- gsub("'", "", packages)
  packages <- gsub(" ", "", packages)

  if (length(packages) == 0) {
    packages <- script[grepl("require", script)]
    packages <- gsub("require\\(", "", packages)
    packages <- gsub("\\)", "", packages)
    packages <- gsub('"', "", packages)
    packages <- gsub("'", "", packages)
    packages <- gsub(" ", "", packages)
  }

  if (length(packages) == 0) {
    packages <- script[grepl("package", script)]
    packages <- gsub("package\\(", "", packages)
    packages <- gsub("\\)", "", packages)
    packages <- gsub('"', "", packages)
    packages <- gsub("'", "", packages)
    packages <- gsub(" ", "", packages)
  }

  output <- capture.output({
    for (pkg in packages) {
      tryCatch(
        {
          utils::install.packages(pkg, dependencies = TRUE)
        },
        error = function(e) {
          print(e)
        }
      )
      tryCatch(
        {
          BiocManager::install(pkg, dependencies = TRUE, update = TRUE, ask = FALSE)
        },
        error = function(e) {
          print(e)
        }
      )
    }
  })
  return(output)
}

#' Drop rows with NA values based on a threshold
#' @description Drop rows with NA values based on a threshold
#' @param data data frame
#' @param percent_allowed_missing numeric, threshold for allowed missing values
#' @return data frame with rows dropped
#' @importFrom dplyr rowwise ungroup filter across everything
drop_na_rows <- function(data, percent_allowed_missing) {
  if (percent_allowed_missing < 0 | percent_allowed_missing > 1) {
    stop("percent_allowed_missing must be between 0 and 1")
  }
  max_na_count <- percent_allowed_missing * ncol(data)
  data %>%
    dplyr::rowwise() %>%
    dplyr::filter(sum(is.na(c_across(everything()))) <= max_na_count) %>%
    dplyr::ungroup()
}

#' Drop columns with NA values based on a threshold
#'
#' @description Drop columns with NA values based on a threshold
#' @param data data frame
#' @param percent_allowed_missing numeric, threshold for allowed missing values
#' @return data frame with columns dropped
#' @importFrom dplyr select where
drop_na_cols <- function(data, percent_allowed_missing) {
  if (percent_allowed_missing < 0 | percent_allowed_missing > 1) {
    stop("percent_allowed_missing must be between 0 and 1")
  }
  max_na_count <- percent_allowed_missing * nrow(data)
  data %>% dplyr::select(where(~ sum(is.na(.)) <= max_na_count))
}

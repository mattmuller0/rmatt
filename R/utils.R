#' @title Utility functions
#' @description Utility functions for general tasks.
#' @name utility_functions
#' @docType package
#' @importFrom utils readLines install.packages
#' @importFrom BiocManager install
#' @importFrom jsonlite write_json
#' @importFrom parallel detectCores
#' @importFrom dplyr rowwise filter ungroup select where
#' @importFrom tidyr drop_na

# LOAD FUNCTIONS
# space reserved for sourcing in functions

###########################################################################
#
#                                 CODE
#
###########################################################################

#' Install all required packages from a script
#'
#' @param script character, path to script
#' @return none
#' @export
install_packages_from_script <- function(script) {
  script <- utils::readLines(script)
  
  packages <- script[grepl('library', script)]
  packages <- gsub('library\\(', '', packages)
  packages <- gsub('suppressMessages\\(', '', packages)
  packages <- gsub('\\)', '', packages)
  packages <- gsub('"', '', packages)
  packages <- gsub("'", '', packages)
  packages <- gsub(' ', '', packages)

  if (length(packages) == 0) {
    packages <- script[grepl('require', script)]
    packages <- gsub('require\\(', '', packages)
    packages <- gsub('\\)', '', packages)
    packages <- gsub('"', '', packages)
    packages <- gsub("'", '', packages)
    packages <- gsub(' ', '', packages)
  }

  if (length(packages) == 0) {
    packages <- script[grepl('package', script)]
    packages <- gsub('package\\(', '', packages)
    packages <- gsub('\\)', '', packages)
    packages <- gsub('"', '', packages)
    packages <- gsub("'", '', packages)
    packages <- gsub(' ', '', packages)
  }

  output <- capture.output({
    for (pkg in packages) {
      tryCatch({
        utils::install.packages(pkg, dependencies = TRUE)
      }, error = function(e) {print(e)})
      tryCatch({
        BiocManager::install(pkg, dependencies = TRUE, update = TRUE, ask = FALSE)
      }, error = function(e) {print(e)})
    }
  })
  return(output)
}

#' Set up a directory for output
#'
#' @param outpath path to output directory
#' @param overwrite logical, whether to overwrite the directory if it already exists
#' @return outpath path to output directory
#' @export
setup_output_dir <- function(outpath, overwrite = FALSE) {
  if (dir.exists(outpath)) {
    if (overwrite) {
      dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
    } else {
      stop("Output directory already exists. Set overwrite = TRUE to overwrite.")
    }
  } else {
    dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
  }
  return(outpath)
}

#' Set the seed for reproducibility
#'
#' @param seed integer, seed to set
#' @return seed integer, seed that was set
#' @export
set_seed <- function(seed) {
  set.seed(seed)
  return(seed)
}

#' Time a function
#'
#' @param func function to time
#' @return time time it took to run the function
#' @export
time_function <- function(func) {
  start_time <- Sys.time()
  func
  end_time <- Sys.time()
  time <- end_time - start_time
  return(time)
}

#' Detect the number of cores on a machine
#'
#' @return cores integer, number of cores on the machine
#' @export
detect_cores <- function() {
  cores <- parallel::detectCores()
  return(cores)
}

#' Detect the number of threads on a machine
#'
#' @return threads integer, number of threads on the machine
#' @export
detect_threads <- function() {
  threads <- parallel::detectCores(logical = FALSE)
  return(threads)
}

#' Generate a JSON file from a list
#'
#' @param list list, list to convert to JSON
#' @param filename character, name of the file to write
#' @return none
#' @export
list_to_json <- function(list, filename) {
  jsonlite::write_json(list, filename)
}

#' Search a vector for a string
#'
#' @param vector vector, vector to search
#' @param string character, string to search for
#' @param ... additional arguments passed to grep
#' @return named vector of indices
#' @export
search_vector <- function(vector, string, ...) {
  idx <- grep(string, vector, ...)
  vector <- vector[idx]
  names(idx) <- vector
  return(idx)
}

#' Drop rows with NA values based on a threshold
#'
#' @param data data frame
#' @param percent_allowed_missing numeric, threshold for allowed missing values
#' @return data frame with rows dropped
#' @export
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
#' @param data data frame
#' @param percent_allowed_missing numeric, threshold for allowed missing values
#' @return data frame with columns dropped
#' @export
drop_na_cols <- function(data, percent_allowed_missing) {
  if (percent_allowed_missing < 0 | percent_allowed_missing > 1) {
    stop("percent_allowed_missing must be between 0 and 1")
  }
  max_na_count <- percent_allowed_missing * nrow(data)
  data %>% dplyr::select(where(~ sum(is.na(.)) <= max_na_count))
}

#' @title Function to check input data
#' @description Function to check input data
#' @param data Data frame containing the data to be checked
#' @param required_columns Vector of required column names
#' @return Logical value indicating whether the data is valid
check_input_data <- function(data, required_columns) {
  # Check if all required columns are present in the data
  missing_columns <- setdiff(required_columns, colnames(data))
  if (length(missing_columns) > 0) {
    stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
  }

  # Check if the data is a data frame
  if (!is.data.frame(data)) {
    stop("Input data must be a data frame.")
  }

  return(TRUE)
}

#' @title Function to calculate eGFR
#' @description Function to calculate estimated glomerular filtration rate (eGFR) according to the CKD-EPI (2021) equation.
#' @param age Age of the individual
#' @param sex vector of sex
#' @param creatinine Serum creatinine level (mg/dL)
#' @return Estimated glomerular filtration rate (eGFR) in mL/min/1.73m²
#' @export
get_eGFR <- function(
  age, sex, creatinine, 
  sex_levels = c("Male", "Female"),
  creatinine_units = c("mg/dL", "umol/L")
) {
  # Validate sex values
  if (!all(sex %in% sex_levels)) {
    stop(paste("Invalid sex values. Must be one of:", paste(sex_levels, collapse = ", ")))
  }
  
  # Validate and standardize creatinine units
  creatinine_units <- match.arg(creatinine_units)
  
  # Convert creatinine to mg/dL if in umol/L
  if (creatinine_units == "umol/L") {
    creatinine <- creatinine * 0.0113
  }
  
  # Calculate eGFR using sex_levels parameter
  ifelse(
    sex == sex_levels[1],
    142 * pmin(creatinine, 1, na.rm = TRUE)**(-0.302) * pmax(creatinine, 1, na.rm = TRUE)**(-1.200) * 0.9938 ** age,
    142 * pmin(creatinine, 1, na.rm = TRUE)**(-0.241) * pmax(creatinine, 1, na.rm = TRUE)**(-1.200) * 0.9938 ** age * 1.012
  )
}

# ======================== Pooled Cohort Functions ========================

# Tests for survival analysis functions
library(testthat)
library(rmatt)
library(survival)
library(ggsurvfit)

# Create mock survival data for testing
create_mock_survival_data <- function(n = 100) {
    set.seed(420)
    data.frame(
        # Basic patient characteristics
        patient_id = 1:n,
        age = rnorm(n, mean = 60, sd = 10),
        sex = factor(sample(c("M", "F"), n, replace = TRUE)),
        group = factor(sample(c("Control", "Treatment"), n, replace = TRUE)),
        multigroup = sample(c("A", "B", "C"), n, replace = TRUE),

        # Survival times and censoring
        time_to_death = rexp(n, rate = 1 / 1000),
        censor_death = sample(c(0, 1), n, replace = TRUE, prob = c(0.3, 0.7)),
        time_to_progression = rexp(n, rate = 1 / 500),
        censor_progression = sample(c(0, 1), n, replace = TRUE, prob = c(0.2, 0.8)),

        # Risk factors
        risk_level = factor(sample(c("Low", "Medium", "High"), n, replace = TRUE)),
        comorbidity = factor(sample(c("Yes", "No"), n, replace = TRUE))
    )
}

# Test log_rank_test function
test_that("log_rank_test produces correct output", {
    # Create mock data
    test_data <- create_mock_survival_data()

    # Test basic functionality
    result <- log_rank_test(
        data = test_data,
        comparisons = "group",
        censors = c("censor_death", "censor_progression"),
        time_prefix = "time_to_",
        censor_prefix = "censor_"
    )

    # Check output structure
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 2)
    expect_equal(ncol(result), 1)
    expect_true(all(result >= 0 & result <= 1)) # p-values should be between 0 and 1

    # Test error handling
    expect_error(log_rank_test(
        data = test_data,
        comparisons = "nonexistent",
        censors = "censor_death"
    ))
})

# Test survival_analysis function
test_that("survival_analysis produces expected output for binary group input", {
    # Create mock data
    test_data <- create_mock_survival_data()
    temp_dir <- tempdir()

    # Test basic functionality
    result <- survival_analysis(
        df = test_data,
        condition = "group",
        censors = c("censor_death", "censor_progression"),
        outdir = temp_dir,
        time_prefix = "time_to_",
        censor_prefix = "censor_"
    )

    # Check output structure
    expect_type(result, "list")
    expect_true(all(c("survival_plots", "hazard_ratios") %in% names(result)))
    expect_true(is.data.frame(result$hazard_ratios))

    # Check if files were created
    expect_true(file.exists(file.path(temp_dir, "hazard_ratios.csv")))
    expect_true(dir.exists(file.path(temp_dir, "survival_plots")))

    # Test error handling
    expect_error(survival_analysis(
        df = test_data,
        condition = "nonexistent",
        censors = c("censor_death"),
        outdir = temp_dir
    ))
})


# Test survival_analysis function
test_that("survival_analysis produces expected output for multigroup input", {
    # Create mock data
    test_data <- create_mock_survival_data()
    temp_dir <- tempdir()

    # Test basic functionality
    result <- survival_analysis(
        df = test_data,
        condition = "multigroup",
        censors = c("censor_death", "censor_progression"),
        outdir = temp_dir,
        time_prefix = "time_to_",
        censor_prefix = "censor_"
    )

    # Check output structure
    expect_type(result, "list")
    expect_true(all(c("survival_plots", "hazard_ratios") %in% names(result)))
    expect_true(is.data.frame(result$hazard_ratios))

    # Check if files were created
    expect_true(file.exists(file.path(temp_dir, "hazard_ratios.csv")))
    expect_true(dir.exists(file.path(temp_dir, "survival_plots")))

    # Test error handling
    expect_error(survival_analysis(
        df = test_data,
        condition = "nonexistent",
        censors = c("censor_death"),
        outdir = temp_dir
    ))
})


# Test hazard_ratios_table function
test_that("hazard_ratios_table handles various input scenarios", {
    # Create mock data
    test_data <- create_mock_survival_data()

    # Test basic functionality
    result <- hazard_ratios_table(
        df = test_data,
        condition = "group",
        censors = c("censor_death", "censor_progression"),
        time_prefix = "time_to_",
        censor_prefix = "censor_"
    )

    # Check output structure
    expect_s3_class(result, "tbl_df") # Check if result is a tibble
    expect_true(all(c("hazard.ratio", "ci.lower", "ci.upper", "p.value") %in% colnames(result)))

    # Test with controls
    result_with_controls <- hazard_ratios_table(
        df = test_data,
        condition = "group",
        censors = c("censor_death"),
        controls = c("age", "sex"),
        time_prefix = "time_to_",
        censor_prefix = "censor_"
    )
    expect_s3_class(result_with_controls, "tbl_df") # Check if result is a tibble

    # Test per_sd option with continuous variable
    result_per_sd <- hazard_ratios_table(
        df = test_data,
        condition = "age",
        censors = c("censor_death"),
        per_sd = TRUE,
        time_prefix = "time_to_",
        censor_prefix = "censor_"
    )
    expect_s3_class(result_per_sd, "tbl_df") # Check if result is a tibble

    # Test one-vs-rest (OVR) analysis
    result_ovr <- hazard_ratios_table(
        df = test_data,
        condition = "risk_level",
        censors = c("censor_death"),
        ovr = TRUE,
        time_prefix = "time_to_",
        censor_prefix = "censor_"
    )
    expect_s3_class(result_ovr, "tbl_df") # Check if result is a tibble
})

# # Test filtered_hazard_ratio_table function
# test_that("filtered_hazard_ratio_table handles filtering correctly", {
#   # Create mock data
#   test_data <- create_mock_survival_data()

#   # Test basic functionality
#   result <- filtered_hazard_ratio_table(
#     data = test_data,
#     condition = "group",
#     risks = c("risk_level", "comorbidity"),
#     censors = c("censor_death", "censor_progression"),
#     time_prefix = "time_to_",
#     censor_prefix = "censor_"
#   )

#   # Check output structure
#   expect_true(is.data.frame(result))
#   expect_true(all(c("x", "y", "n") %in% colnames(result)))

#   # Test error handling for non-character risks
#   expect_error(filtered_hazard_ratio_table(
#     data = test_data,
#     condition = "group",
#     risks = list(1, 2), # Invalid risks type
#     censors = c("censor_death")
#   ))

#   # Test with verbose option
#   result_verbose <- filtered_hazard_ratio_table(
#     data = test_data,
#     condition = "group",
#     risks = c("risk_level"),
#     censors = c("censor_death"),
#     verbose = TRUE,
#     time_prefix = "time_to_",
#     censor_prefix = "censor_"
#   )
#   expect_true(is.data.frame(result_verbose))
# })

# # Test edge cases and error handling
# test_that("functions handle edge cases appropriately", {
#   # Create mock data
#   test_data <- create_mock_survival_data()

#   # Test with missing data
#   test_data$time_to_death[1] <- NA
#   test_data$censor_death[2] <- NA
#   expect_message(hazard_ratios_table(
#     df = test_data,
#     condition = "group",
#     censors = c("censor_death"),
#     time_prefix = "time_to_",
#     censor_prefix = "censor_"
#   ))

#   # Test with single observation per group
#   small_data <- test_data[1:2, ]
#   expect_error(log_rank_test(
#     data = small_data,
#     comparisons = "group",
#     censors = c("censor_death")
#   ), NA)

#   # Test with all censored data
#   all_censored_data <- test_data
#   all_censored_data$censor_death <- 0
#   expect_warning(survival_analysis(
#     df = all_censored_data,
#     condition = "group",
#     censors = c("censor_death"),
#     outdir = tempdir()
#   ), NA)
# })

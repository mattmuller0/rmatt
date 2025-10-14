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
        outpath = temp_dir,
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
        outpath = temp_dir
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
        outpath = temp_dir,
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
        outpath = temp_dir
    ))
})

# Test hazard_ratios_table function
test_that("hazard_ratios_table works correctly with basic input", {
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
    expect_true(is.data.frame(result))
    expect_true(all(c("censor", "condition", "term", "hazard.ratio", 
                     "ci.lower", "ci.upper", "p.value") %in% colnames(result)))
    expect_equal(nrow(result), 2)  # One row per censor
    expect_true(all(result$hazard.ratio > 0))  # HRs should be positive
})

test_that("hazard_ratios_table handles controls correctly", {
    test_data <- create_mock_survival_data()
    
    result <- hazard_ratios_table(
        df = test_data,
        condition = "group",
        censors = "censor_death",
        controls = c("age", "sex"),
        time_prefix = "time_to_",
        censor_prefix = "censor_"
    )
    
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 1)  # One row for the condition
})

test_that("hazard_ratios_table handles per_sd option correctly", {
    test_data <- create_mock_survival_data()
    
    # Test with continuous variable
    result <- hazard_ratios_table(
        df = test_data,
        condition = "age",
        censors = "censor_death",
        per_sd = TRUE,
        time_prefix = "time_to_",
        censor_prefix = "censor_"
    )
    
    expect_true(is.data.frame(result))
    
    # Should error with categorical variable
    expect_error(
        hazard_ratios_table(
            df = test_data,
            condition = "group",
            censors = "censor_death",
            per_sd = TRUE
        ),
        "Condition must be numeric for per_sd"
    )
})

test_that("hazard_ratios_table handles tibble input correctly", {
    test_data <- create_mock_survival_data()
    test_data <- tidyr::as_tibble(test_data)
    
    result <- hazard_ratios_table(
        df = test_data,
        condition = "multigroup",
        censors = "censor_death",
        ovr = TRUE,
        time_prefix = "time_to_",
        censor_prefix = "censor_"
    )
    
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 3)  # One row per level of multigroup
})

test_that("hazard_ratios_table handles one-vs-rest correctly", {
    test_data <- create_mock_survival_data()
    
    result <- hazard_ratios_table(
        df = test_data,
        condition = "multigroup",
        censors = "censor_death",
        ovr = TRUE,
        time_prefix = "time_to_",
        censor_prefix = "censor_"
    )
    
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 3)  # One row per level of multigroup
})

test_that("hazard_ratios_table handles NA values appropriately", {
    test_data <- create_mock_survival_data()
    test_data$age[1:10] <- NA
    
    expect_warning(
        result <- hazard_ratios_table(
            df = test_data,
            condition = "age",
            censors = "censor_death",
            time_prefix = "time_to_",
            censor_prefix = "censor_"
        ),
        "Removing 10 rows with NA values"
    )
    
    expect_true(is.data.frame(result))
})

test_that("hazard_ratios_table validates input columns", {
    test_data <- create_mock_survival_data()
    
    expect_error(
        hazard_ratios_table(
            df = test_data,
            condition = "nonexistent",
            censors = "censor_death"
        )
    )
    
    expect_error(
        hazard_ratios_table(
            df = test_data,
            condition = "group",
            censors = "nonexistent"
        )
    )
})

test_that("hazard_ratios_table handles subgroup analyses correctly", {
    test_data <- create_mock_survival_data()
    test_data$subgroup <- factor(sample(c("X", "Y"), nrow(test_data), replace = TRUE))
    
    result <- hazard_ratios_table(
        df = test_data,
        condition = "group",
        censors = "censor_death",
        subgroups = "subgroup",
        time_prefix = "time_to_",
        censor_prefix = "censor_"
    )
    
    expect_true(is.data.frame(result))
    expect_true(all(c("subgroup_var", "subgroup_val") %in% colnames(result)))
})

# Tests for multiple conditions functionality
test_that("log_rank_test works with multiple conditions", {
    test_data <- create_mock_survival_data()
    
    # Test multiple conditions (always univariate)
    result <- log_rank_test(
        data = test_data,
        comparisons = c("group", "sex"),
        censors = c("censor_death", "censor_progression")
    )
    
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 2)  # 2 censors
    expect_equal(ncol(result), 2)  # 2 conditions
    expect_true(all(colnames(result) %in% c("group", "sex")))
})

test_that("survival_analysis works with multiple conditions", {
    test_data <- create_mock_survival_data()
    temp_dir <- tempdir()
    
    # Test multiple conditions (always univariate)
    result <- survival_analysis(
        df = test_data,
        condition = c("group", "sex"),
        censors = c("censor_death", "censor_progression"),
        outpath = temp_dir
    )
    
    expect_type(result, "list")
    expect_true(all(c("survival_plots", "hazard_ratios") %in% names(result)))
    expect_true(is.data.frame(result$hazard_ratios))
    expect_true("condition_var" %in% colnames(result$hazard_ratios))
    
    # Check that separate plots were created for each condition
    expected_plots <- c("group_censor_death", "group_censor_progression", 
                       "sex_censor_death", "sex_censor_progression")
    expect_true(all(expected_plots %in% names(result$survival_plots)))
})

test_that("hazard_ratios_table works with multiple conditions", {
    test_data <- create_mock_survival_data()
    
    # Test multiple conditions (always univariate)
    result <- hazard_ratios_table(
        df = test_data,
        condition = c("group", "sex"),
        censors = c("censor_death", "censor_progression")
    )
    
    expect_true(is.data.frame(result))
    expect_true("condition" %in% colnames(result))
    unique_conditions <- unique(result$condition)
    expect_true(all(c("group", "sex") %in% unique_conditions))
})

test_that("multiple conditions functionality maintains backward compatibility", {
    test_data <- create_mock_survival_data()
    temp_dir <- tempdir()
    
    # Test that single condition (as vector) works the same as before
    result_single_vec <- survival_analysis(
        df = test_data,
        condition = c("group"),  # Single condition as vector
        censors = c("censor_death"),
        outpath = temp_dir
    )
    
    result_single_char <- survival_analysis(
        df = test_data,
        condition = "group",  # Single condition as character
        censors = c("censor_death"),
        outpath = temp_dir
    )
    
    # Both should produce similar structure (allowing for some differences in naming)
    expect_type(result_single_vec, "list")
    expect_type(result_single_char, "list")
    expect_true(all(c("survival_plots", "hazard_ratios") %in% names(result_single_vec)))
    expect_true(all(c("survival_plots", "hazard_ratios") %in% names(result_single_char)))
})
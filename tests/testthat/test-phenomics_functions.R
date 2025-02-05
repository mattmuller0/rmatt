# tests/testthat/test-phenomics.R

test_that("make_composite_coding works correctly", {
  codes <- list(
    ICD10 = c("I21", "I22"),
    OPCS = c("K45", "K75")
  )
  
  result <- make_composite_coding(codes, "MI")
  
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 3)
  expect_equal(unique(result$key), "MI")
  expect_equal(sort(unique(result$field)), sort(names(codes)))
})

test_that("make_phenotypes processes data correctly", {
  hesin <- data.frame(
    "Participant ID" = rep(1:3, each = 2),
    field = rep(c("ICD10", "OPCS"), 3),
    value = c("I21", "K45", "I22", "K75", "I23", "K46"),
    stringsAsFactors = FALSE
  )
  
  encoding <- data.frame(
    key = "MI",
    field = c("ICD10", "ICD10", "OPCS"),
    value = c("I21", "I22", "K45"),
    stringsAsFactors = FALSE
  )
  
  result <- make_phenotypes(hesin, encoding)
  
  expect_true("phenotype" %in% colnames(result))
  expect_true("n_codes" %in% colnames(result))
  expect_equal(unique(result$phenotype), "MI")
})

test_that("functions handle edge cases", {
  # Empty input
  expect_error(make_composite_coding(list(), "test"))
  
  # Missing columns
  hesin_bad <- data.frame(ID = 1, bad = "ICD10", value = "I21")
  encoding_good <- data.frame(key = "MI", field = "ICD10", value = "I21")
  expect_error(make_phenotypes(hesin_bad, encoding_good))
  
  # No matches
  hesin_good <- data.frame(
    ID = 1,
    field = "ICD10",
    value = "I21"
  )
  encoding_nomatch <- data.frame(
    key = "MI",
    field = "ICD10",
    value = "I99"
  )
  result <- make_phenotypes(hesin_good, encoding_nomatch)
  expect_equal(nrow(result), 0)
})
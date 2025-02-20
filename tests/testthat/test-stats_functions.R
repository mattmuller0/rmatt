# tests/testthat/test-stats-functions.R

test_that("stats_table creates correct table", {
  data <- data.frame(
    group = factor(rep(c("A", "B"), each = 50)),
    age = c(rnorm(50, 60, 10), rnorm(50, 65, 10)),
    sex = factor(sample(c("M", "F"), 100, replace = TRUE)),
    bmi = c(rnorm(50, 25, 3), rnorm(50, 27, 3))
  )
  
  result <- stats_table(data, "group", c("age", "sex", "bmi"))
  print(result)
  
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 6)  # Labels + level + Overall + 2 groups + p-value
  expect_true(all(c("Overall", "A", "B", "p") %in% colnames(result)))
})

test_that("hypergeometric_scoring calculates correctly", {
  # Mock GSEA object
  MockGSEA <- setClass("MockGSEA", 
    slots = c(
      gene = "character",
      geneSets = "list",
      universe = "character",
      result = "data.frame"
    )
  )
  
  gse <- new("MockGSEA",
    gene = c("gene1", "gene2", "gene3"),
    geneSets = list(
      set1 = c("gene1", "gene4", "gene5"),
      set2 = c("gene2", "gene6", "gene7")
    ),
    universe = c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene7"),
    result = data.frame(ID = c("set1", "set2"))
  )
  
  result <- hypergeometric_scoring(gse, method = "fisher")
  
  expect_true("odds_ratio" %in% colnames(result@result))
  expect_equal(length(result@result$odds_ratio), 2)
})

test_that("correlation_matrix produces valid correlations", {
  data <- data.frame(
    var1 = rnorm(100),
    var2 = rnorm(100),
    var3 = rnorm(100)
  )
  
  result <- correlation_matrix(
    data, 
    vars1 = c("var1", "var2"),
    vars2 = c("var2", "var3")
  )
  
  expect_true(all(c("cor_matr", "pvalue_matr", "cor_matr_filtr") %in% names(result)))
  expect_true(all(abs(result$cor_matr) <= 1))
  expect_true(all(result$pvalue_matr >= 0 & result$pvalue_matr <= 1))
})

test_that("correlation_matrix handles missing values", {
  data <- data.frame(
    var1 = c(1:9, NA),
    var2 = c(NA, 2:10),
    var3 = 1:10
  )
  
  result <- correlation_matrix(
    data,
    vars1 = c("var1", "var2"),
    vars2 = c("var2", "var3"),
    use = "pairwise.complete.obs"
  )
  
  expect_false(any(is.na(result$cor_matr)))
  expect_false(any(is.na(result$pvalue_matr)))
})

test_that("hypergeometric_scoring validates method", {
  MockGSEA <- setClass("MockGSEA", 
    slots = c(
      gene = "character",
      geneSets = "list",
      universe = "character",
      result = "data.frame"
    )
  )
  
  gse <- new("MockGSEA",
    gene = c("gene1"),
    geneSets = list(set1 = c("gene1")),
    universe = c("gene1"),
    result = data.frame(ID = "set1")
  )
  
  expect_error(hypergeometric_scoring(gse, method = "invalid"))
  expect_error(hypergeometric_scoring(gse, method = NULL))
})
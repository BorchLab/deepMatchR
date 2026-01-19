# tests/testthat/test-toSerology.R

# --- Single Allele Conversion Tests ---

test_that("toSerology converts Class I alleles correctly", {
  # A locus
  expect_equal(toSerology("A*01:01"), "A1")
  expect_equal(toSerology("A*02:01"), "A2")
  expect_equal(toSerology("A*03:01"), "A3")

  # B locus
  expect_equal(toSerology("B*07:02"), "B7")
  expect_equal(toSerology("B*08:01"), "B8")

  # C locus (note: returns "Cw" prefix)
  result <- toSerology("C*01:02")
  expect_true(grepl("^Cw", result))
})

test_that("toSerology converts Class II alleles correctly", {
  # DRB1 locus
  expect_equal(toSerology("DRB1*03:01"), "DR17")
  expect_equal(toSerology("DRB1*04:01"), "DR4")

  # DQB1 locus
  result <- toSerology("DQB1*02:01")
  expect_true(grepl("^DQ", result))
})

test_that("toSerology handles DRB3/4/5 loci",
{
  # DRB3 -> DR52
  result <- toSerology("DRB3*01:01")
  expect_true(grepl("^DR", result))

  # DRB4 -> DR53
  result <- toSerology("DRB4*01:01")
  expect_true(grepl("^DR", result))

  # DRB5 -> DR51
  result <- toSerology("DRB5*01:01")
  expect_true(grepl("^DR", result))
})

# --- Vector Input Tests ---

test_that("toSerology handles vector input", {
  alleles <- c("A*01:01", "B*07:02", "DRB1*03:01")
  result <- toSerology(alleles)

  expect_length(result, 3)
  expect_equal(result[1], "A1")
  expect_equal(result[2], "B7")
  expect_equal(result[3], "DR17")
})

test_that("toSerology handles NA/empty values in vectors", {
  alleles <- c("A*01:01", NA, "", "B*07:02")
  result <- toSerology(alleles)

  expect_length(result, 4)
  expect_equal(result[1], "A1")
  expect_true(is.na(result[2]))
  expect_true(is.na(result[3]))
  expect_equal(result[4], "B7")
})

# --- hla_genotype Input Tests ---

test_that("toSerology works with hla_genotype - serology return", {
  df <- data.frame(
    A_1 = "A*01:01", A_2 = "A*02:01",
    B_1 = "B*07:02", B_2 = "B*08:01",
    stringsAsFactors = FALSE
  )
  geno <- hlaGeno(df)

  result <- toSerology(geno, return = "serology")

  expect_true(is.character(result))
  expect_true("A1" %in% result)
  expect_true("A2" %in% result)
  expect_true("B7" %in% result)
  expect_true("B8" %in% result)
})

test_that("toSerology works with hla_genotype - genotype return", {
  df <- data.frame(
    A_1 = "A*01:01", A_2 = "A*02:01",
    stringsAsFactors = FALSE
  )
  geno <- hlaGeno(df)

  result <- toSerology(geno, return = "genotype")

  expect_s3_class(result, "hla_genotype")
  expect_true("A_ser_1" %in% names(result$data))
  expect_true("A_ser_2" %in% names(result$data))
  expect_equal(result$data$A_ser_1, "A1")
  expect_equal(result$data$A_ser_2, "A2")
})

test_that("toSerology works with hla_genotype - data.frame return", {
  df <- data.frame(
    A_1 = "A*01:01", A_2 = "A*02:01",
    stringsAsFactors = FALSE
  )
  geno <- hlaGeno(df)

  result <- toSerology(geno, return = "data.frame")

  expect_s3_class(result, "data.frame")
  expect_true(all(c("allele", "locus", "serology", "serology_full") %in% names(result)))
})

# --- data.frame Return Tests ---

test_that("toSerology data.frame return has correct structure", {
  result <- toSerology(c("A*01:01", "B*07:02"), return = "data.frame")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(c("allele", "locus", "allele_2f", "serology", "serology_full") %in% names(result)))
})

# --- Split Resolution Tests ---

test_that("toSerology respects resolve_splits parameter", {
  # DRB1*15:01 should be DR15 (split of DR2) when resolve_splits = TRUE
  with_split <- toSerology("DRB1*15:01", resolve_splits = TRUE)
  without_split <- toSerology("DRB1*15:01", resolve_splits = FALSE)

  # Both should start with DR

  expect_true(grepl("^DR", with_split))
  expect_true(grepl("^DR", without_split))
})

# --- Multi-field Allele Tests ---

test_that("toSerology handles multi-field alleles", {
  # Three-field
  expect_equal(toSerology("A*01:01:01"), "A1")

  # Four-field
  expect_equal(toSerology("A*01:01:01:01"), "A1")
})

# --- Error Handling Tests ---

test_that("toSerology validates input types", {
  expect_error(toSerology(123), "character")
  expect_error(toSerology(list(a = "A*01:01")), "character")
})

test_that("toSerology respects na_action parameter", {
  # NA action - should return NA silently
  result <- toSerology("INVALID*99:99", na_action = "NA")
  expect_true(is.na(result))

  # warn action - should return NA with warning
  expect_warning(
    toSerology("INVALID*99:99", na_action = "warn"),
    "No serology mapping"
  )

  # error action - should stop
  expect_error(
    toSerology("INVALID*99:99", na_action = "error"),
    "No serology mapping"
  )
})

test_that("toSerology handles malformed alleles gracefully", {
  # Missing asterisk
  result <- toSerology("A0101", na_action = "NA")
  expect_true(is.na(result))

  # Wrong format
  result <- toSerology("not_an_allele", na_action = "NA")
  expect_true(is.na(result))
})

# --- Edge Cases ---

test_that("toSerology handles empty input", {
  result <- toSerology(character(0))
  expect_length(result, 0)
})

test_that("toSerology handles single NA input", {
  result <- toSerology(NA_character_)
  expect_true(is.na(result))
})

test_that("toSerology handles whitespace in alleles", {
  expect_equal(toSerology("  A*01:01  "), "A1")
})

# --- updateWmdaData Tests ---

test_that("clearWmdaCache works without error", {
  # Should not error even if no cache exists
  expect_silent(clearWmdaCache(verbose = FALSE))
})

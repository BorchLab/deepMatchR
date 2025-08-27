# tests/testthat/test-quantifyMismatch.R

test_that("basic mismatch counting works", {
  seq1 <- "YFAMYGEKVAHTHVDTLYVRYHY"
  seq2 <- "YFDMYGEKVAHTHVDTLYVRFHY"
  expect_equal(quantifyMismatch(seq1, seq2), 2L)
})

test_that("case-insensitivity and equal length checks", {
  expect_error(quantifyMismatch("ACD", "ACDE"), "same length")
  expect_error(quantifyMismatch(1, "ACD"), "character strings")
  
  # Same result regardless of case
  expect_equal(quantifyMismatch("acde", "AcDe"), 0L)
})

test_that("filters: polarity-only vs charge-only", {
  # Construct sequences with three mismatches:
  # 1) A -> D (nonpolar -> negative)   => charge_change=TRUE,  polarity_change=TRUE
  # 2) A -> S (nonpolar -> polar)      => charge_change=FALSE, polarity_change=TRUE
  # 3) S -> T (polar -> polar)         => charge_change=FALSE, polarity_change=FALSE (still mismatch)
  s1 <- "AAAS"
  s2 <- "ADST"
  
  # Raw mismatches: 3
  expect_equal(quantifyMismatch(s1, s2), 3L)
  
  # Only charge-changing: counts #1
  expect_equal(quantifyMismatch(s1, s2, filter_charge = TRUE), 1L)
  
  # Only polarity-changing: counts #1 and #2
  expect_equal(quantifyMismatch(s1, s2, filter_polarity = TRUE), 2L)
  
  # Require BOTH charge and polarity change: counts #1
  expect_equal(
    quantifyMismatch(s1, s2, filter_charge = TRUE, filter_polarity = TRUE),
    1L
  )
  
  # Only mismatches that do NOT change polarity: counts #3
  expect_equal(quantifyMismatch(s1, s2, filter_polarity = FALSE), 1L)
})

test_that("charge change without polarity change is possible (charged->charged)", {
  # E (neg, polar) -> D (neg, polar): mismatch but charge_change=FALSE, polarity_change=FALSE
  expect_equal(
    quantifyMismatch("E", "D", filter_charge = TRUE),
    0L
  )
  
  # C (polar, neutral) -> E (polar, negative): charge_change=TRUE, polarity_change=FALSE (both polar)
  expect_equal(
    quantifyMismatch("C", "E", filter_charge = TRUE, filter_polarity = FALSE),
    1L
  )
})

test_that("na_action behavior with unknowns (X, *)", {
  s1 <- "ACDX"
  s2 <- "ACDY"
  
  # error: complains about unknown residues
  expect_error(quantifyMismatch(s1, s2, na_action = "error"), "Unknown/unsupported")
  
  # exclude: unknown-involving positions not counted when filters apply
  # Here only the X vs Y position is a mismatch, but it's unknown on s1
  expect_equal(quantifyMismatch(s1, s2, na_action = "exclude"), 1L) # no filters => still 1
  # With a filter active, the unknown position is excluded (NA) from counting
  expect_equal(
    quantifyMismatch(s1, s2, na_action = "exclude", filter_polarity = TRUE),
    0L
  )
  
  # count: mismatches vs unknowns are counted; property deltas set to NA internally
  expect_equal(quantifyMismatch(s1, s2, na_action = "count"), 1L)
})

test_that("return types and columns are correct", {
  s1 <- "ACDE"
  s2 <- "ACDF"
  
  # Count
  expect_type(quantifyMismatch(s1, s2, return = "count"), "integer")
  
  # data.frame
  df <- quantifyMismatch(s1, s2, return = "detail")
  expect_s3_class(df, "data.frame")
  expect_true(all(c(
    "position","ref","alt","is_mismatch",
    "charge_ref","charge_alt","charge_change",
    "polarity_ref","polarity_alt","polarity_change","counted"
  ) %in% names(df)))
  
})

test_that("identical sequences return zero", {
  s <- "MSTNPKPQR"
  expect_equal(quantifyMismatch(s, s), 0L)
})


context("Testing Eplet Mismatch Quantification")

test_that("quantifyEpletMismatch correctly counts mismatches", {
  data(deepMatchR_eplets)
  allele1 <- "A*01:110" 
  allele2 <- "A*02:636" 
  expect_equal(quantifyEpletMismatch(allele1, allele2), 12)
})

test_that("quantifyEpletMismatch handles no mismatches", {
  expect_equal(quantifyEpletMismatch("A*01:110", "A*01:110"), 0)
})

test_that("quantifyEpletMismatch handles alleles not in the database", {
  # One allele not in db
  expect_equal(quantifyEpletMismatch("A*01:110", "A*99:99"), 9) 

  # Both alleles not in db
  expect_equal(quantifyEpletMismatch("A*98:98", "A*99:99"), 0)
})

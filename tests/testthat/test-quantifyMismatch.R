context("Testing Mismatch Quantification")

test_that("quantifyMismatch correctly counts mismatches", {
  seq1 <- "ABCDE"
  seq2 <- "ABXDE"
  expect_equal(quantifyMismatch(seq1, seq2), 1)

  seq3 <- "ABCDE"
  seq4 <- "XXXXX"
  expect_equal(quantifyMismatch(seq3, seq4), 5)

  seq5 <- "ABCDE"
  seq6 <- "ABCDE"
  expect_equal(quantifyMismatch(seq5, seq6), 0)
})

test_that("quantifyMismatch handles errors correctly", {
  expect_error(quantifyMismatch("ABC", "ABCD"), "Input sequences must be of the same length.")
  expect_error(quantifyMismatch(123, "ABC"), "Input sequences must be character strings.")
})

test_that("getAlleleSequence retrieves a known allele", {
  skip_if_not_installed("immReferent")
  # Skip test if IMGT is not available (e.g., no internet)
  if (!immReferent::is_imgt_available()) {
    skip("IMGT website not available.")
  }

  # A*01:01 is a very common allele, it should exist.
  seq <- getAlleleSequence("A*01:01")
  expect_type(seq, "character")
  expect_gt(nchar(seq), 0) # Expect that the sequence is not empty
})

test_that("getAlleleSequence throws an error for a non-existent allele", {
  skip_if_not_installed("immReferent")
  if (!immReferent::is_imgt_available()) {
    skip("IMGT website not available.")
  }

  expect_error(
    getAlleleSequence("A*99:99"),
    "Allele 'A\\*99:99' not found in the IMGT/HLA database."
  )
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

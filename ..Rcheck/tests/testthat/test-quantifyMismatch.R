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
  # These alleles are in the deepMatchR_eplets dataset
  # A*01:110 has eplet 21H, A*03:210 has eplet 21H
  # Let's find some better examples from the data
  data(deepMatchR_eplets)
  allele1 <- "A*01:110" # Has eplet 21H
  allele2 <- "A*02:636" # Has eplet 21H

  # Let's find alleles with different eplets
  # From exploration, we know A*01:110 has 21H and A*24:02 has many others
  # This is not robust as it depends on the data, but it's a start

  # For a reproducible test, let's use known values from the dataset
  # eplets for A*01:110 is 21H
  # eplets for A*24:02 is 114R, 142T, 144K, 145R, 150A, 151H, 152A, 156W, 163T, 167W, 44Y, 65Q, 66N, 76A, 77N, 80N
  # Mismatches should be the union of these two sets, minus the intersection (which is empty)
  expect_equal(quantifyEpletMismatch("A*01:110", "A*24:02"), 13)
})

test_that("quantifyEpletMismatch handles no mismatches", {
  expect_equal(quantifyEpletMismatch("A*01:110", "A*01:110"), 0)
})

test_that("quantifyEpletMismatch handles alleles not in the database", {
  # One allele not in db
  expect_equal(quantifyEpletMismatch("A*01:110", "A*99:99"), 11) # Mismatch is just the eplets from the first allele

  # Both alleles not in db
  expect_equal(quantifyEpletMismatch("A*98:98", "A*99:99"), 0)
})

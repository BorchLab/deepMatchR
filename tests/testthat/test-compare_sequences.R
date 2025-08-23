context("Testing compare_hla_sequences function")

test_that("compare_hla_sequences handles identical sequences", {
  seq1 <- "YFAMYGEKVAHTHVDTLYVRYHY"
  seq2 <- "YFAMYGEKVAHTHVDTLYVRYHY"
  result <- compare_hla_sequences(seq1, seq2)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
})

test_that("compare_hla_sequences stops with unequal length sequences", {
  seq1 <- "YFAMYGEKVAHTHVDTLYVRYHY"
  seq2 <- "YFAMYGEKVAHTHVDTLYVRYH"
  expect_error(compare_hla_sequences(seq1, seq2), "Input sequences must be of the same length.")
})

test_that("compare_hla_sequences identifies polymorphisms correctly", {
  seq1 <- "YFAMYGEKVAHTHVDTLYVRYHY"
  seq2 <- "YFDMYGEKVAHTHVDTLYVRFHY"
  result <- compare_hla_sequences(seq1, seq2)

  expect_equal(nrow(result), 2)

  # Check first polymorphism
  expect_equal(result$Position[1], 3)
  expect_equal(result$AA_Seq1[1], "A")
  expect_equal(result$AA_Seq2[1], "D")
  expect_equal(result$Polarity_Change[1], "Nonpolar to Polar")
  expect_equal(result$Charge_Change[1], "Uncharged to Negative")

  # Check second polymorphism
  expect_equal(result$Position[2], 21)
  expect_equal(result$AA_Seq1[2], "Y")
  expect_equal(result$AA_Seq2[2], "F")
  expect_equal(result$Polarity_Change[2], "Polar to Nonpolar")
  expect_equal(result$Charge_Change[2], "No change")
})

test_that("compare_hla_sequences handles unknown amino acids", {
  seq1 <- "XFAMYGEKVAHTHVDTLYVRYHY"
  seq2 <- "YFAMYGEKVAHTHVDTLYVRYHY"
  expect_warning(
    result <- compare_hla_sequences(seq1, seq2),
    "Unknown amino acid at position 1. Skipping property analysis for this position."
  )
  expect_equal(nrow(result), 1)
  expect_equal(result$Position[1], 1)
  expect_equal(result$AA_Seq1[1], "X")
  expect_equal(result$AA_Seq2[1], "Y")
  expect_equal(result$Polarity_Change[1], "Unknown")
  expect_equal(result$Charge_Change[1], "Unknown")

  seq1_2 <- "YFAMYGEKVAHTHVDTLYVRYHY"
  seq2_2 <- "YFAMYGEKVZHTHVDTLYVRYHY"
    expect_warning(
    result2 <- compare_hla_sequences(seq1_2, seq2_2),
    "Unknown amino acid at position 10. Skipping property analysis for this position."
  )
  expect_equal(nrow(result2), 1)
  expect_equal(result2$Position[1], 10)
  expect_equal(result2$AA_Seq1[1], "A")
  expect_equal(result2$AA_Seq2[1], "Z")
  expect_equal(result2$Polarity_Change[1], "Unknown")
  expect_equal(result2$Charge_Change[1], "Unknown")
})

test_that("compare_hla_sequences handles non-character inputs", {
    expect_error(compare_hla_sequences(123, "abc"), "Input sequences must be character strings.")
    expect_error(compare_hla_sequences("abc", 123), "Input sequences must be character strings.")
    expect_error(compare_hla_sequences(NULL, "abc"), "Input sequences must be character strings.")
})

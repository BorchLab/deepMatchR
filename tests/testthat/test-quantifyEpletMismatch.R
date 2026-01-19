
# tests/testthat/test-quantifyEpletMismatch.R

test_that("quantifyEpletMismatch returns 0 for identical alleles", {
  result <- quantifyEpletMismatch("A*01:01", "A*01:01")
  expect_equal(result, 0L)
})

test_that("quantifyEpletMismatch returns integer", {
  result <- quantifyEpletMismatch("A*01:01", "A*02:01")
  expect_type(result, "integer")
  expect_true(result >= 0)
})

test_that("quantifyEpletMismatch validates input types", {
  expect_error(quantifyEpletMismatch(123, "A*01:01"), "character")
  expect_error(quantifyEpletMismatch("A*01:01", 456), "character")
})
test_that("quantifyEpletMismatch validates input length", {
  expect_error(
    quantifyEpletMismatch(c("A*01:01", "A*02:01"), "A*03:01"),
    "length-1"
  )
  expect_error(
    quantifyEpletMismatch("A*01:01", c("A*02:01", "A*03:01")),
    "length-1"
  )
})

test_that("quantifyEpletMismatch is symmetric", {
  result_ab <- quantifyEpletMismatch("A*01:01", "A*02:01")
  result_ba <- quantifyEpletMismatch("A*02:01", "A*01:01")
  expect_equal(result_ab, result_ba)
})

test_that("quantifyEpletMismatch respects evidence_level filter", {
  # Broader filter should give >= result than narrower filter
  all_evidence <- quantifyEpletMismatch("A*01:01", "A*02:01",
                                         evidence_level = c("A1", "A2", "B", "D"))
  a1_a2_only <- quantifyEpletMismatch("A*01:01", "A*02:01",
                                       evidence_level = c("A1", "A2"))
  a1_only <- quantifyEpletMismatch("A*01:01", "A*02:01",
                                    evidence_level = "A1")

  expect_true(a1_only <= a1_a2_only)
  expect_true(a1_a2_only <= all_evidence)
})

test_that("quantifyEpletMismatch works with NULL evidence_level", {
  result <- quantifyEpletMismatch("A*01:01", "A*02:01", evidence_level = NULL)
  expect_type(result, "integer")
  expect_true(result >= 0)
})

test_that("quantifyEpletMismatch works across different loci", {
  # A vs A
  result_a <- quantifyEpletMismatch("A*01:01", "A*02:01")
  expect_true(result_a >= 0)

  # B vs B
  result_b <- quantifyEpletMismatch("B*07:02", "B*08:01")
  expect_true(result_b >= 0)
})
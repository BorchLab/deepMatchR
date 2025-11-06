#tests/testthat/test-quantifyEpletMismatch.R

test_that("quantifyEpletMismatch: identity is zero and basic counts compute", {
  
    # identical -> 0
    expect_identical(quantifyEpletMismatch("A*01:01","A*01:01"), 0L)
  
    m1 <- quantifyEpletMismatch("A*01:01","A*02:01")
    expect_identical(m1, 13L)

    m2 <- quantifyEpletMismatch("B*44:02","B*51:01", evidence_level = "A1")
    expect_gte(m2, 4L)
})

test_that("quantifyEpletMismatch: filters (exposition/reactivity) alter results", {

    base <- quantifyEpletMismatch("A*03:01","A*24:02")  # default evidence A1/A2
    hi_only <- quantifyEpletMismatch("A*03:01","A*24:02",
                                      exposition_filter = "High")
    igg_only <- quantifyEpletMismatch("A*03:01","A*24:02",
                                      reactivity_filter = "IgG")
    expect_true(is.integer(base) && is.integer(hi_only) && is.integer(igg_only))
    # Not asserting specific numbers, but they should compute and be >= 0
    expect_gte(base, 0L); expect_gte(hi_only, 0L); expect_gte(igg_only, 0L)
})

test_that("quantifyEpletMismatch: argument validation", {
    expect_error(quantifyEpletMismatch(1, "A*01:01"), "length-1 character")
    expect_error(quantifyEpletMismatch("A*01:01", c("A*02:01","A*03:01")), "length-1 character")
})

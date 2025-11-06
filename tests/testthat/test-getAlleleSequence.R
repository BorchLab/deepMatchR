# tests/testthat/test-getAlleleSequence.R

mock_seq_db <- list(
  "A*01:01" = "YFAMYGEKVAHTHVDTLYVRYHY", # 23 AA
  "A*02:01" = "YFDMYGEKVAHTHVDTLYVRFHY", # 23 AA (Mismatch at 3 (A->D), 21 (Y->F))
  "A*03:01" = "YFAMYGEKVAHTHVDTLYVRYXX", # 23 AA (Mismatch at 22, 23 (unknown))
  "B*07:02" = "GSHSMRYFYTAMSRPGRGEPRFI", # 23 AA
  "B*08:01" = "GSHSMRYFYTAMSRPGRGEPRFI", # 23 AA (Identical to B*07:02)
  "B*44:02" = "GSHSMRYFYTAMSREGRGEPRFI", # 23 AA (Mismatch at 14 (P->E))
  "B*51:01" = "GSHSMRYFYTAMSREGRGEPRFX", # 23 AA (Mismatch at 14 (P->E), 23 (unknown))
  "C*01:02" = "AAAAAAAAAA",
  "C*02:02" = "AAAAAAAAAA",
  "C*03:04" = "AABBAAAAAA" # Mismatch at 3, 4
)

test_that("getAlleleSequence and batchGetSequences work", {
  
  # Mock immReferent::getIMGT
  mock_imgt <- function(gene, type, suppressMessages = TRUE) {
    if (!identical(type, "PROT")) {
      stop("mock_imgt only supports type = 'PROT' in this test")
    }
    # Return a named character vector: names = alleles, values = sequences
    c(
      "A*01:01" = mock_seq_db[["A*01:01"]],
      "A*02:01" = mock_seq_db[["A*02:01"]],
      "B*07:02" = mock_seq_db[["B*07:02"]]
    )
  }
  
  testthat::with_mocked_bindings(
    # Bindings that live in your test's environment (if any)
    # .getAlleleSequenceBase = .getAlleleSequenceBase,  # likely not needed
    
    # Bindings that live in the immReferent namespace
    getIMGT = mock_imgt,
    .package = "immReferent",
    {
      expect_equal(getAlleleSequence("A*01:01", type = "PROT", use_cache = FALSE),
                   mock_seq_db[["A*01:01"]])
      
      alleles <- c("A*01:01", "B*07:02", "A*01:01")
      expected_list <- list(
        "A*01:01" = mock_seq_db[["A*01:01"]],
        "B*07:02" = mock_seq_db[["B*07:02"]]
      )
      
      batch_res <- batchGetSequences(alleles, type = "PROT", n_cores = 1, use_cache = FALSE)
      expect_equal(length(batch_res), 2)
      expect_equal(batch_res, expected_list)
      
      expect_error(
        getAlleleSequence("Z*99:99", type = "PROT", use_cache = FALSE),
        "Allele 'Z\\*99:99' not found"
      )
    })
})
  
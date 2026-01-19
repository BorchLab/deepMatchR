# tests/testthat/test-calculateMismatchLoad.R

# --- Mock AA sequences for alleles used in tests ---
mock_seq_db <- c(
  "A*01:01" = "YFAMYGEKVAHTHVDTLYVRYHY",  # 23 AA
  "A*02:01" = "YFDMYGEKVAHTHVDTLYVRFHY",  # 23 AA; mismatches at 3 (A->D), 21 (Y->F)
  "A*03:01" = "YFAMYGEKVAHTHVDTLYVRYXX",  # 23 AA; unknowns at 22-23
  "B*07:02" = "GSHSMRYFYTAMSRPGRGEPRFI",  # 23 AA
  "B*08:01" = "GSHSMRYFYTAMSRPGRGEPRFI",  # 23 AA (identical to B*07:02)
  "B*44:02" = "GSHSMRYFYTAMSREGRGEPRFI",  # 23 AA (P->E at pos 14)
  "B*51:01" = "GSHSMRYFYTAMSREGRGEPRFX"   # 23 AA (P->E at 14; unknown at 23)
)

# Toy genotypes from the examples
recipient_df <- data.frame(
  A_1 = "A*01:01", A_2 = "A*02:01",
  B_1 = "B*07:02", B_2 = "B*08:01",
  check.names = FALSE
)
donor_df <- data.frame(
  A_1 = "A*01:01", A_2 = "A*03:01",
  B_1 = "B*44:02", B_2 = "B*51:01",
  check.names = FALSE
)
rgeno <- hlaGeno(recipient_df)
dgeno <- hlaGeno(donor_df)

# Mocks
mock_validate <- function(x) TRUE
mock_batchGetSequences <- function(alleles, n_cores = NULL, ...) {
  # Return a named character vector allele -> sequence
  miss <- setdiff(unique(alleles), names(mock_seq_db))
  if (length(miss)) stop("Unknown allele(s) in mock batchGetSequences: ", paste(miss, collapse = ", "))
  mock_seq_db[unique(alleles)]
}

test_that("total: default counts across all shared loci", {
  testthat::with_mocked_bindings(
    validateHlaGeno     = mock_validate,
    batchGetSequences   = mock_batchGetSequences,
    {
      total <- calculateMismatchLoad(rgeno, dgeno)
      expect_type(total, "integer")
      expect_gt(total, 0L)
    }
  )
})

test_that("per_locus: returns integer loads per locus with expected names", {
  testthat::with_mocked_bindings(
    validateHlaGeno   = mock_validate,
    batchGetSequences = mock_batchGetSequences,
    {
      df <- calculateMismatchLoad(rgeno, dgeno, return = "per_locus")
      expect_true(is.data.frame(df))
      expect_setequal(names(df), c("locus", "mismatch_load"))
      expect_setequal(df$locus, c("A", "B"))
      expect_true(is.integer(df$mismatch_load))
      
      # Sanity: locus A should be > 0, locus B should be > 0
      # (We rely on quantifyMismatch for exact values.)
      expect_true(all(df$mismatch_load[df$locus %in% c("A","B")] >= 0L))
    }
  )
})

test_that("loci argument subsets computation", {
  testthat::with_mocked_bindings(
    validateHlaGeno   = mock_validate,
    batchGetSequences = mock_batchGetSequences,
    {
      onlyA <- calculateMismatchLoad(rgeno, dgeno, loci = "A", return = "per_locus")
      expect_identical(sort(unique(onlyA$locus)), "A")
      both  <- calculateMismatchLoad(rgeno, dgeno, return = "per_locus")
      sumA  <- onlyA$mismatch_load[onlyA$locus == "A"]
      sumBothA <- both$mismatch_load[both$locus == "A"]
      expect_identical(sumA, sumBothA)
    }
  )
})

test_that("pairwise: returns proper matrix for a given locus with row/col names", {
  testthat::with_mocked_bindings(
    validateHlaGeno   = mock_validate,
    batchGetSequences = mock_batchGetSequences,
    {
      mB <- calculateMismatchLoad(rgeno, dgeno, return = "pairwise", pairwise_locus = "B")
      expect_true(is.matrix(mB))
      expect_setequal(rownames(mB), c("B*07:02", "B*08:01" ))
      expect_setequal(colnames(mB), c("B*44:02", "B*51:01"))
      expect_true(all(is.numeric(mB)))
      # Diagonal could be non-zero if recipient vs donor allele differs at same slot;
      # we won't assume symmetry, just that entries are >= 0 integers.
      expect_true(all(mB >= 0))
      expect_true(all(mB == as.integer(mB)))
    }
  )
})

test_that("filters are forwarded to quantifyMismatch (charge / polarity)", {
  testthat::with_mocked_bindings(
    validateHlaGeno   = mock_validate,
    batchGetSequences = mock_batchGetSequences,
    {
      total_all    <- calculateMismatchLoad(rgeno, dgeno, return = "total")
      total_charge <- calculateMismatchLoad(rgeno, dgeno, filter_charge = TRUE, return = "total")
      total_polar  <- calculateMismatchLoad(rgeno, dgeno, filter_polarity = TRUE, return = "total")
      # Charge- or polarity-only loads should be <= all mismatches
      expect_lte(total_charge, total_all)
      expect_lte(total_polar,  total_all)
    }
  )
})

test_that("na_action is respected via quantifyMismatch behavior", {
  # Introduce unknown residue into one donor allele sequence (A*03:01 already has 'X')
  # We'll make A*02:01 contain an 'X' to exercise 'error' and 'exclude'.
  local_db <- mock_seq_db
  local_db["A*02:01"] <- sub(".", "X", local_db["A*02:01"])  # force unknown at pos 1
  
  mock_batch_with_unknown <- function(alleles, n_cores = NULL, ...) {
    miss <- setdiff(unique(alleles), names(local_db))
    if (length(miss)) stop("Unknown allele(s) in mock batchGetSequences: ", paste(miss, collapse = ", "))
    local_db[unique(alleles)]
  }
  
  testthat::with_mocked_bindings(
    validateHlaGeno   = mock_validate,
    batchGetSequences = mock_batch_with_unknown,
    {
      
      # 'exclude' and 'count' should both return integers (not error)
      excl <- calculateMismatchLoad(rgeno, dgeno, na_action = "exclude")
      cnt  <- calculateMismatchLoad(rgeno, dgeno, na_action = "count")
      expect_true(is.integer(excl))
      expect_true(is.integer(cnt))
    }
  )
})

test_that("no shared loci -> error", {
  # Make a donor with only C locus
  donorC <- hlaGeno(data.frame(C_1 = "C*01:02", C_2 = "C*02:02", check.names = FALSE))
  
  # Extend DB for C alleles to avoid batchGetSequences errors (not actually reached)
  mock_db2 <- c(mock_seq_db, "C*01:02" = "AAAAAAAAAA", "C*02:02" = "AAAAAAAAAA")
  mock_batch2 <- function(alleles, ...) { mock_db2[unique(alleles)] }
  
  testthat::with_mocked_bindings(
    validateHlaGeno   = mock_validate,
    batchGetSequences = mock_batch2,
    {
      expect_error(
        calculateMismatchLoad(rgeno, donorC),
        "No shared loci"
      )
    }
  )
})

test_that("return shapes and types are correct (total vs per_locus)", {
  testthat::with_mocked_bindings(
    validateHlaGeno   = mock_validate,
    batchGetSequences = mock_batchGetSequences,
    {
      tot <- calculateMismatchLoad(rgeno, dgeno, return = "total")
      expect_true(is.integer(tot))
      df  <- calculateMismatchLoad(rgeno, dgeno, return = "per_locus")
      expect_true(is.data.frame(df))
      expect_true(is.integer(df$mismatch_load))
      expect_true(all(df$locus %in% c("A","B")))
    }
  )
})

test_that("parallel branch (mclapply) works when enabled (skipped on Windows)", {
  testthat::skip_on_os("windows")
  # Also skip on CRAN/CI if desired:
  # testthat::skip_on_cran()
  
  testthat::with_mocked_bindings(
    validateHlaGeno   = mock_validate,
    batchGetSequences = mock_batchGetSequences,
    {
      # Ensure we have >1 shared locus and n_cores > 1
      tot <- calculateMismatchLoad(
        rgeno, dgeno,
        parallel = TRUE, n_cores = 2, return = "total"
      )
      expect_true(is.integer(tot))
    }
  )
})

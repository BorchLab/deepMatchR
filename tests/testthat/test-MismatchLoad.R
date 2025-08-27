# tests/testthat/test-calculateMismatchLoad.R

# Minimal toy allele->sequence map (equal length strings).
.local_seq_map <- c(
  "A*01:01" = "AAAAAA",
  "A*02:01" = "AAAIAA", # 1 mismatch vs A*01:01 (pos4 I vs A)
  "A*03:01" = "AAADAA", # 1 mismatch vs A*01:01 (pos4 D vs A)  (charge change + polarity change)
  "B*07:02" = "CCCCCC",
  "B*08:01" = "CCCACC", # 1 mismatch vs B*07:02 (pos4 A vs C) (nonpolar->nonpolar)
  "B*44:02" = "CCCECC", # 1 mismatch vs B*07:02 (pos4 E vs C) (neutral/polar C -> negative/polar E): charge change TRUE, polarity change FALSE
  "B*51:01" = "CCCCRC"  # 1 mismatch vs B*07:02 (pos5 R vs C) (neutral/polar C -> positive/polar R): charge change TRUE, polarity change FALSE
)

# Example genotypes
make_example_genos <- function() {
  recipient <- data.frame(
    A_1 = "A*01:01", A_2 = "A*02:01",
    B_1 = "B*07:02", B_2 = "B*08:01",
    stringsAsFactors = FALSE
  )
  donor <- data.frame(
    A_1 = "A*01:01", A_2 = "A*03:01",
    B_1 = "B*44:02", B_2 = "B*51:01",
    stringsAsFactors = FALSE
  )
  list(r = hlaGeno(recipient), d = hlaGeno(donor))
}

test_that("total load equals manual sum over loci", {
  gens <- make_example_genos()
  
  # Manual computation with mocked getAlleleSequence
  testthat::with_mocked_bindings(
    {
      # Build all pairwise counts per locus using quantifyMismatch("count")
      get <- function(a) .local_seq_map[[a]]
      # A locus pairs:
      A_r <- c("A*01:01","A*02:01"); A_d <- c("A*01:01","A*03:01")
      total_A <- 0L
      for (ra in A_r) for (da in A_d) {
        if (ra == da) next
        total_A <- total_A + quantifyMismatch(get(ra), get(da), return = "count")
      }
      # B locus pairs:
      B_r <- c("B*07:02","B*08:01"); B_d <- c("B*44:02","B*51:01")
      total_B <- 0L
      for (ra in B_r) for (da in B_d) {
        if (ra == da) next
        total_B <- total_B + quantifyMismatch(get(ra), get(da), return = "count")
      }
      manual_total <- total_A + total_B
      
      # Function result
      fun_total <- calculateMismatchLoad(gens$r, gens$d)
      
      expect_equal(fun_total, manual_total)
    },
    getAlleleSequence = function(allele) .local_seq_map[[allele]]
  )
})

test_that("per_locus returns a sensible data.frame and sums to total", {
  gens <- make_example_genos()
  
  testthat::with_mocked_bindings(
    {
      per <- calculateMismatchLoad(gens$r, gens$d, return = "per_locus")
      expect_s3_class(per, "data.frame")
      expect_true(all(c("locus", "mismatch_load") %in% names(per)))
      expect_true(all(per$locus %in% c("A","B")))
      expect_equal(calculateMismatchLoad(gens$r, gens$d),
                   sum(per$mismatch_load))
    },
    getAlleleSequence = function(allele) .local_seq_map[[allele]]
  )
})



test_that("loci subsetting restricts computation", {
  gens <- make_example_genos()
  
  testthat::with_mocked_bindings(
    {
      only_A <- calculateMismatchLoad(gens$r, gens$d, loci = "A")
      only_B <- calculateMismatchLoad(gens$r, gens$d, loci = "B")
      both   <- calculateMismatchLoad(gens$r, gens$d)
      
      expect_equal(only_A + only_B, both)
    },
    getAlleleSequence = function(allele) .local_seq_map[[allele]]
  )
})

test_that("filters are forwarded to quantifyMismatch", {
  gens <- make_example_genos()
  
  testthat::with_mocked_bindings(
    {
      # Because our B-locus donor alleles include charged residues (E, R),
      # charge-only counting should reduce or equal the unfiltered total.
      total_all    <- calculateMismatchLoad(gens$r, gens$d)
      total_charge <- calculateMismatchLoad(gens$r, gens$d, filter_charge = TRUE)
      total_pol    <- calculateMismatchLoad(gens$r, gens$d, filter_polarity = TRUE)
      
      expect_true(total_charge <= total_all)
      expect_true(total_pol    <= total_all)
      
      # Per-locus with filters should sum to the filtered total
      per_charge <- calculateMismatchLoad(gens$r, gens$d, return = "per_locus", filter_charge = TRUE)
      expect_equal(sum(per_charge$mismatch_load), total_charge)
    },
    getAlleleSequence = function(allele) .local_seq_map[[allele]]
  )
})

test_that("pairwise requires a valid locus and shared loci must exist", {
  gens <- make_example_genos()
  
  testthat::with_mocked_bindings(
    {
      # Missing pairwise_locus argument
      expect_error(
        calculateMismatchLoad(gens$r, gens$d, return = "pairwise"),
        "pairwise_locus"
      )
      
      # Non-shared locus
      expect_error(
        calculateMismatchLoad(gens$r, gens$d, return = "pairwise", pairwise_locus = "C"),
        "not in shared_loci"
      )
    },
    getAlleleSequence = function(allele) .local_seq_map[[allele]]
  )
  
  # No shared loci case
  r_only_A <- hlaGeno(data.frame(A_1 = "A*01:01", stringsAsFactors = FALSE))
  d_only_B <- hlaGeno(data.frame(B_1 = "B*07:02", stringsAsFactors = FALSE))
  expect_error(
    calculateMismatchLoad(r_only_A, d_only_B),
    "No shared loci"
  )
})

test_that("handles NA/empty allele cells gracefully", {
  r <- hlaGeno(data.frame(A_1 = "A*01:01", A_2 = NA_character_, stringsAsFactors = FALSE))
  d <- hlaGeno(data.frame(A_1 = "A*03:01", A_2 = "",            stringsAsFactors = FALSE))
  
  testthat::with_mocked_bindings(
    {
      # Should compute using present strings only; not error on NA/""
      val <- calculateMismatchLoad(r, d, loci = "A")
      expect_true(is.integer(val) || is.numeric(val))
    },
    getAlleleSequence = function(allele) .local_seq_map[[allele]]
  )
})

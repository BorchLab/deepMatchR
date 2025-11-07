# tests/testthat/test-compareLoadMetrics.R

test_that("compareLoadMetrics merges sources and computes derived columns", {
  # Mock return tables
  fake_mismatch <- data.frame(
    locus = c("A","B","C"),
    mismatch_load = c(2, 0, 5)
  )
  fake_eplet <- data.frame(
    locus = c("A","B","C"),
    eplet_load = c(4, 1, 2)
  )
  fake_pep <- data.frame(
    locus = c("A","B","C"),
    binding_peptides   = c(10, 0, 20),
    binding_percentage = c(25, 0, 50),
    stringsAsFactors = FALSE
  )
  
  # capture passthrough args
  seen <- new.env(parent = emptyenv())
  seen$loci <- NULL
  seen$mhc_class <- NULL
  
  with_mocked_bindings(
    calculateMismatchLoad = function(recipient_geno, donor_geno, loci, return) {
      seen$loci <- loci
      expect_identical(return, "per_locus")
      fake_mismatch
    },
    calculateEpletLoad = function(recipient_geno, donor_geno, loci, return) {
      expect_identical(loci, seen$loci)
      expect_identical(return, "per_locus")
      fake_eplet
    },
    calculatePeptideBindingLoad = function(recipient_geno, donor_geno, loci, mhc_class, return) {
      seen$mhc_class <- mhc_class
      expect_identical(loci, seen$loci)
      expect_identical(return, "per_locus")
      fake_pep
    },
    {
      out <- compareLoadMetrics(
        recipient_geno = list(),
        donor_geno = list(),
        loci = c("A","B","C"),
        mhc_class = "I"
      )
      
      # structure
      expect_true(all(c("locus","mismatch_load","eplet_load",
                        "binding_peptides","binding_percentage",
                        "mismatch_eplet_ratio","binding_per_mismatch") %in% names(out)))
      expect_equal(nrow(out), 3L)
      
      # passthrough confirmation
      expect_identical(seen$loci, c("A","B","C"))
      expect_identical(seen$mhc_class, "I")
      
      # derived columns:
      # mismatch_eplet_ratio = mismatch / (eplet + 1)
      expect_equal(out$mismatch_eplet_ratio[out$locus == "A"], 2/(4+1))
      expect_equal(out$mismatch_eplet_ratio[out$locus == "B"], 0/(1+1))
      expect_equal(out$mismatch_eplet_ratio[out$locus == "C"], 5/(2+1))
      
      # binding_per_mismatch = binding_peptides / (mismatch + 1)
      expect_equal(out$binding_per_mismatch[out$locus == "A"], 10/(2+1))
      expect_equal(out$binding_per_mismatch[out$locus == "B"], 0/(0+1))
      expect_equal(out$binding_per_mismatch[out$locus == "C"], 20/(5+1))
    }
  )
})

test_that("compareLoadMetrics works if loci is NULL (pass-through)", {
  with_mocked_bindings(
    calculateMismatchLoad = function(recipient_geno, donor_geno, loci, return) {
      expect_null(loci)
      data.frame(locus = "A", mismatch_load = 1)
    },
    calculateEpletLoad = function(recipient_geno, donor_geno, loci, return) {
      expect_null(loci)
      data.frame(locus = "A", eplet_load = 2)
    },
    calculatePeptideBindingLoad = function(recipient_geno, donor_geno, loci, mhc_class, return) {
      expect_null(loci)
      data.frame(locus = "A", binding_peptides = 3, binding_percentage = 50)
    },
    {
      out <- compareLoadMetrics(list(), list(), loci = NULL, mhc_class = "II")
      expect_equal(out$locus, "A")
      expect_equal(out$mismatch_eplet_ratio, 1/(2+1))
      expect_equal(out$binding_per_mismatch, 3/(1+1))
    }
  )
})

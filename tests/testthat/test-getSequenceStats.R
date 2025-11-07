# tests/testthat/test-getSequenceStats.R


test_that("protein stats are computed correctly", {
  seqs <- c(
    A = "ACCPPWYV",     # 2 Cys, 2 Pro, hydrophobic = A,V,I,L,M,F,Y,W -> here A,V,W,Y = 4/8
    B = "GGGG----**",   # gaps/star ignored
    C = "stnqyr"        # lower-case and mixed letters (will be uppercased)
  )
  res <- getSequenceStats(seqs, type = "PROT")
  expect_equal(res$allele, names(seqs))
  expect_equal(res$n_cysteines[res$allele == "A"], 2L)
  expect_equal(res$n_prolines[res$allele == "A"], 2L)
  expect_equal(res$length[res$allele == "A"], as.vector(nchar(seqs["A"])))
  # hydrophobic ratio for A: A,V,W,Y = 4/8 = 0.5
  expect_equal(res$hydrophobic_ratio[res$allele == "A"], 0.5, tolerance = 1e-12)
  
  # Gaps and '*' ignored -> no invalid counted there; but raw length still includes them
  # invalid count for B should be 0 after ignoring ignore_chars
  expect_equal(res$n_invalid[res$allele == "B"], 0L)
  
  # Lower-case handled
  expect_true(res$frac_aromatic[res$allele == "C"] >= 0 && res$frac_aromatic[res$allele == "C"] <= 1)
  expect_false(any(is.na(res$kmer_entropy))) # should compute something
})

test_that("nucleotide content and stops work as intended", {
  seqs <- c(
    X = "ATGC",       # GC=0.5, AT=0.5
    Y = "ATGTAA",     # has stop (TAA) in frame 0
    Z = "NNNNACGT"    # ambiguous Ns should be excluded from GC/AT denom
  )
  res <- getSequenceStats(seqs, type = "NUC")
  expect_equal(res$gc_content[res$allele == "X"], 0.5, tolerance = 1e-12)
  expect_equal(res$at_content[res$allele == "X"], 0.5, tolerance = 1e-12)
  expect_true(res$has_stop_codon_anyframe[res$allele == "Y"])
  # For Z, only ACGT count toward denom; "ACGT" -> GC=0.5 AT=0.5
  expect_equal(res$gc_content[res$allele == "Z"], 0.5, tolerance = 1e-12)
  expect_equal(res$at_content[res$allele == "Z"], 0.5, tolerance = 1e-12)
})

test_that("reference comparisons behave correctly", {
  seqs <- c(S1 = "ACDEFG", S2 = "ACREFG", S3 = "ACDE")
  # Use S1 as reference by name
  res <- getSequenceStats(seqs, type = "PROT", ref = "S1")
  # S1 vs S1 = identity 1, Hamming 0
  expect_equal(res$identity_to_ref[res$allele == "S1"], 1.0, tolerance = 1e-12)
  expect_equal(res$hamming_to_ref[res$allele == "S1"], 0)
  
  # S2 differs at position 3 (X vs D) -> identity 5/6 ≈ 0.8333; Hamming 1
  expect_equal(res$identity_to_ref[res$allele == "S2"], 5/6, tolerance = 1e-12)
  expect_equal(res$hamming_to_ref[res$allele == "S2"], 1)
  
  # S3 shorter; identity computed over shortest (4) -> first four all match
  expect_equal(res$identity_to_ref[res$allele == "S3"], 1.0, tolerance = 1e-12)
  expect_true(!is.na(res$hamming_to_ref[res$allele == "S3"]))
})

test_that("raw reference sequence is accepted", {
  seqs <- c(A = "AAAAAA", B = "AAAATA")
  res <- getSequenceStats(seqs, type = "NUC", ref = "AAAAAA")
  expect_equal(res$identity_to_ref[res$allele == "A"], 1.0, tolerance = 1e-12)
  # B differs at position 5 -> identity 5/6
  expect_equal(res$identity_to_ref[res$allele == "B"], 5/6, tolerance = 1e-12)
  expect_equal(res$ref_label[1], "<raw_ref>")
})

test_that("pairwise identity matrix is symmetric and sane", {
  seqs <- c(A = "ACGT", B = "ACGA", C = "TCGT")
  out <- getSequenceStats(seqs, type = "NUC", compute_pairs = TRUE)
  M <- out$pairwise_identity
  expect_true(all(diag(M) == 1))
  expect_equal(M["A","B"], M["B","A"])
  # A vs B differ at last pos -> 3/4
  expect_equal(M["A","B"], 3/4, tolerance = 1e-12)
})

test_that("invalid characters are counted and warning can be silenced", {
  seqs <- c(A = "AZZAC", B = "ACGTN-")
  expect_warning(
    res <- getSequenceStats(seqs, type = "NUC"),
    regexp = "Invalid symbols detected"
  )
  # Silenced
  res2 <- getSequenceStats(seqs, type = "NUC", warn_invalid = FALSE)
  expect_s3_class(res2, "data.frame")
  # '-' ignored by default; 'Z' is invalid for NUC and should be counted
  expect_true(res2$n_invalid[res2$allele == "A"] >= 1L)
})

test_that("k-mer entropy returns NA for too-short sequences", {
  seqs <- c(A = "A", B = "AC")
  # With k=3, both too short
  res_p <- getSequenceStats(seqs, type = "PROT", k = 3)
  expect_true(all(is.na(res_p$kmer_entropy)))
  res_n <- getSequenceStats(seqs, type = "NUC", k = 3)
  expect_true(all(is.na(res_n$kmer_entropy)))
})

test_that("type arg is validated", {
  seqs <- c(A = "AC")
  expect_error(getSequenceStats(seqs, type = "XYZ"), "arg")
})

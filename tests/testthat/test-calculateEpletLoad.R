# tests/testthat/test-calculateEpletLoad.R

# Example genotypes for testing
make_eplet_genos <- function() {
  recipient <- data.frame(
    A_1 = "A*01:01", A_2 = "A*02:01",
    B_1 = "B*07:02", B_2 = "B*08:01",
    stringsAsFactors = FALSE
  )
  donor <- data.frame(
    A_1 = "A*03:01", A_2 = "A*24:02",
    B_1 = "B*44:02", B_2 = "B*51:01",
    stringsAsFactors = FALSE
  )
  list(r = hlaGeno(recipient), d = hlaGeno(donor))
}

test_that("calculateEpletLoad returns integer for total", {
  gens <- make_eplet_genos()
  result <- calculateEpletLoad(gens$r, gens$d)

  expect_type(result, "integer")
  expect_true(result >= 0)
})

test_that("calculateEpletLoad per_locus returns correct structure", {
  gens <- make_eplet_genos()
  per_locus <- calculateEpletLoad(gens$r, gens$d, return = "per_locus")

  expect_s3_class(per_locus, "data.frame")
  expect_true(all(c("locus", "eplet_load") %in% names(per_locus)))
  expect_true(all(per_locus$locus %in% c("A", "B")))
})

test_that("calculateEpletLoad per_locus sums to total", {
  gens <- make_eplet_genos()

  total <- calculateEpletLoad(gens$r, gens$d)
  per_locus <- calculateEpletLoad(gens$r, gens$d, return = "per_locus")

  expect_equal(sum(per_locus$eplet_load), total)
})

test_that("calculateEpletLoad pairwise returns matrix", {
  gens <- make_eplet_genos()

  mat <- calculateEpletLoad(gens$r, gens$d, return = "pairwise", pairwise_locus = "A")

  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 2)  # recipient alleles
  expect_equal(ncol(mat), 2)  # donor alleles
  expect_true(all(mat >= 0))
})

test_that("calculateEpletLoad respects evidence_level filter", {
  gens <- make_eplet_genos()

  # Broader filter should give >= result than narrower filter
  all_evidence <- calculateEpletLoad(gens$r, gens$d, evidence_level = c("A1", "A2", "B", "D"))
  a1_a2_only <- calculateEpletLoad(gens$r, gens$d, evidence_level = c("A1", "A2"))
  a1_only <- calculateEpletLoad(gens$r, gens$d, evidence_level = "A1")

  expect_true(a1_only <= a1_a2_only)
  expect_true(a1_a2_only <= all_evidence)
})

test_that("calculateEpletLoad respects loci filter", {
  gens <- make_eplet_genos()

  only_A <- calculateEpletLoad(gens$r, gens$d, loci = "A")
  only_B <- calculateEpletLoad(gens$r, gens$d, loci = "B")
  both <- calculateEpletLoad(gens$r, gens$d)

  expect_equal(only_A + only_B, both)
})

test_that("calculateEpletLoad returns 0 for identical genotypes", {
  r <- hlaGeno(data.frame(A_1 = "A*01:01", A_2 = "A*02:01", stringsAsFactors = FALSE))
  d <- hlaGeno(data.frame(A_1 = "A*01:01", A_2 = "A*02:01", stringsAsFactors = FALSE))

  result <- calculateEpletLoad(r, d)
  expect_equal(result, 0L)
})

test_that("calculateEpletLoad errors on no shared loci", {
  r <- hlaGeno(data.frame(A_1 = "A*01:01", stringsAsFactors = FALSE))
  d <- hlaGeno(data.frame(B_1 = "B*07:02", stringsAsFactors = FALSE))

  expect_error(calculateEpletLoad(r, d), "No shared loci")
})

test_that("calculateEpletLoad validates genotype inputs", {
  gens <- make_eplet_genos()

  expect_error(calculateEpletLoad("not_a_geno", gens$d), "hla_genotype")
  expect_error(calculateEpletLoad(gens$r, "not_a_geno"), "hla_genotype")
})

test_that("calculateEpletLoad pairwise requires valid locus", {
  gens <- make_eplet_genos()

  # Missing pairwise_locus
  expect_error(
    calculateEpletLoad(gens$r, gens$d, return = "pairwise"),
    "pairwise_locus"
  )

  # Non-shared locus
  expect_error(
    calculateEpletLoad(gens$r, gens$d, return = "pairwise", pairwise_locus = "C"),
    "not shared"
  )
})

test_that("calculateEpletLoad handles NA/empty allele cells", {
  r <- hlaGeno(data.frame(A_1 = "A*01:01", A_2 = NA_character_, stringsAsFactors = FALSE))
  d <- hlaGeno(data.frame(A_1 = "A*03:01", A_2 = "", stringsAsFactors = FALSE))

  # Should compute using present strings only
  val <- calculateEpletLoad(r, d, loci = "A")
  expect_true(is.integer(val))
})

test_that("calculateEpletLoad honors loci restriction & filters", {
  rgeno <- hlaGeno(data.frame(A_1="A*01:01", A_2="A*02:01", B_1="B*07:02"))
  dgeno <- hlaGeno(data.frame(A_1="A*03:01", A_2="A*24:02", B_1="B*44:02"))
  
  # Only locus A
  tot_A <- calculateEpletLoad(rgeno, dgeno, loci = "A", return = "total")
  # Only locus B
  tot_B <- calculateEpletLoad(rgeno, dgeno, loci = "B", return = "total")
  expect_true(is.integer(tot_A) && is.integer(tot_B))
  expect_gte(tot_A, 0L); expect_gte(tot_B, 0L)
      
  # Evidence filter reduces/changes counts
  tot_default <- calculateEpletLoad(rgeno, dgeno, return = "total")
  tot_A1_only <- calculateEpletLoad(rgeno, dgeno, evidence_level = "A1", return = "total")
  expect_true(is.integer(tot_default))
  expect_true(is.integer(tot_A1_only))
})

test_that("calculateEpletLoad: input validation & error paths", {
  # No shared loci after filter
  rgeno <- hlaGeno(data.frame(A_1="A*01:01"))
  dgeno <- hlaGeno(data.frame(B_1="B*44:02"))
  
  expect_error(calculateEpletLoad(rgeno, dgeno, loci = "A"), "No shared loci")
  
  # Missing alleles (empty strings) -> error
  rgeno2 <- hlaGeno(data.frame(A_1 = "", A_2 = NA))
  dgeno2 <- hlaGeno(data.frame(A_1 = "A*03:01"))
  expect_error(calculateEpletLoad(rgeno2, dgeno2), "No allele strings")

  
  # Pairwise: bad locus or not shared
  rgeno3 <- hlaGeno(data.frame(A_1="A*01:01", A_2="A*02:01"))
  dgeno3 <- hlaGeno(data.frame(A_1="A*03:01", A_2="A*24:02"))
  expect_error(calculateEpletLoad(rgeno3, dgeno3, return="pairwise", pairwise_locus = 1),
                 "provide pairwise_locus")
  expect_error(calculateEpletLoad(rgeno3, dgeno3, return="pairwise", pairwise_locus = "B"),
                 "not shared")
})

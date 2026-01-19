# tests/testthat/test-calculatePeptideBindingLoad.R

# Example genotypes for testing
make_binding_genos <- function() {
  recipient <- data.frame(
    A_1 = "A*02:01", A_2 = "A*03:01",
    stringsAsFactors = FALSE
  )
  donor <- data.frame(
    A_1 = "A*01:01", A_2 = "A*24:02",
    stringsAsFactors = FALSE
  )
  list(r = hlaGeno(recipient), d = hlaGeno(donor))
}

test_that("calculatePeptideBindingLoad accepts raw peptides", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "NLVPMVATV", "FLKEKGGL")

  result <- calculatePeptideBindingLoad(recipient, peptides)

  expect_type(result, "double")
  expect_true(result >= 0)
})

test_that("calculatePeptideBindingLoad returns numeric for total", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "NLVPMVATV")

  result <- calculatePeptideBindingLoad(recipient, peptides, return = "total")

  expect_type(result, "double")
})

test_that("calculatePeptideBindingLoad returns data.frame for detail", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "NLVPMVATV")

  result <- calculatePeptideBindingLoad(recipient, peptides, return = "detail")

  expect_s3_class(result, "data.frame")
  expect_true(all(c("peptide", "hla_allele", "predicted_ic50", "binding_level", "contribution") %in% names(result)))
})

test_that("calculatePeptideBindingLoad returns data.frame for summary", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", A_2 = "A*03:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "NLVPMVATV")

  result <- calculatePeptideBindingLoad(recipient, peptides, return = "summary")

  expect_s3_class(result, "data.frame")
  expect_true(all(c("hla_allele", "n_peptides", "n_strong", "n_weak", "risk_contribution") %in% names(result)))
  expect_equal(nrow(result), 2)  # Two recipient alleles
})

test_that("calculatePeptideBindingLoad validates recipient input", {
  peptides <- c("GILGFVFTL")

  expect_error(
    calculatePeptideBindingLoad("not_valid", peptides),
    "hla_genotype|allele"
  )
})

test_that("calculatePeptideBindingLoad returns 0 for empty peptides", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- character(0)

  result <- calculatePeptideBindingLoad(recipient, peptides, return = "total")
  expect_equal(result, 0)
})

test_that("calculatePeptideBindingLoad filters peptides by length", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "SHORT", "TOOLONGPEPTIDE")  # 9, 5, 14

  result <- calculatePeptideBindingLoad(recipient, peptides, peptide_length = 9L, return = "detail")

  # Should only include 9-mers
  expect_true(all(nchar(result$peptide) == 9))
})

test_that("calculatePeptideBindingLoad aggregate_method works", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "NLVPMVATV", "FLKEKGGL")

  result_sum <- calculatePeptideBindingLoad(recipient, peptides, aggregate_method = "sum")
  result_max <- calculatePeptideBindingLoad(recipient, peptides, aggregate_method = "max")
  result_mean <- calculatePeptideBindingLoad(recipient, peptides, aggregate_method = "mean")

  expect_true(result_max <= result_sum)
  expect_true(result_mean <= result_sum)
})

test_that("calculatePeptideBindingLoad binding thresholds work", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL")

  # With stricter threshold, should have lower risk
  result_strict <- calculatePeptideBindingLoad(
    recipient, peptides,
    binding_threshold = 100, weak_threshold = 500
  )

  result_loose <- calculatePeptideBindingLoad(
    recipient, peptides,
    binding_threshold = 1000, weak_threshold = 10000
  )

  # Both should be numeric
  expect_type(result_strict, "double")
  expect_type(result_loose, "double")
})

test_that("calculatePeptideBindingLoad classifies binding levels", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL")  # Known A*02:01 binder

  result <- calculatePeptideBindingLoad(
    recipient, peptides,
    return = "detail",
    binding_threshold = 500,
    weak_threshold = 5000
  )

  # binding_level should be one of: strong, weak, non_binder
  expect_true(all(result$binding_level %in% c("strong", "weak", "non_binder")))
})

test_that("calculatePeptideBindingLoad PWM backend works with multiple alleles", {
  recipient <- hlaGeno(data.frame(
    A_1 = "A*02:01", A_2 = "A*03:01",
    B_1 = "B*07:02", B_2 = "B*44:02",
    stringsAsFactors = FALSE
  ))
  peptides <- c("GILGFVFTL", "NLVPMVATV")

  result <- calculatePeptideBindingLoad(recipient, peptides, backend = "pwm", return = "detail")

  # Should have predictions for each peptide x allele combination
  expect_equal(nrow(result), length(peptides) * 4)  # 2 peptides x 4 alleles
})

test_that("calculatePeptideBindingLoad netmhcpan requires backend_path", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL")

  expect_error(
    calculatePeptideBindingLoad(recipient, peptides, backend = "netmhcpan"),
    "backend_path"
  )
})

test_that("calculatePeptideBindingLoad accepts character vector of alleles", {
  alleles <- c("A*02:01", "A*03:01")
  peptides <- c("GILGFVFTL", "NLVPMVATV")

  result <- calculatePeptideBindingLoad(alleles, peptides)

  expect_type(result, "double")
  expect_true(result >= 0)
})

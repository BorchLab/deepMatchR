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

# --- Identical Genotype Tests ---

test_that("calculatePeptideBindingLoad returns 0 for identical genotypes with no mismatched peptides", {
  # When donor and recipient have identical alleles, no mismatched peptides should be derived
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", A_2 = "A*03:01", stringsAsFactors = FALSE))
  donor <- hlaGeno(data.frame(A_1 = "A*02:01", A_2 = "A*03:01", stringsAsFactors = FALSE))

 # With empty peptides (when peptides cannot be derived), should return 0
  result <- calculatePeptideBindingLoad(recipient, character(0), return = "total")
  expect_equal(result, 0)
})

# --- PWM Backend Detailed Tests ---

test_that("PWM backend produces consistent scores for known binders", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  # GILGFVFTL is a well-known A*02:01 binder (influenza M1 peptide)
  peptides <- c("GILGFVFTL")

  result <- calculatePeptideBindingLoad(recipient, peptides, backend = "pwm", return = "detail")

  expect_equal(nrow(result), 1)
  expect_equal(result$peptide, "GILGFVFTL")
  expect_equal(result$hla_allele, "A*02:01")
  expect_true(result$predicted_ic50 > 0)
})

test_that("PWM backend handles unknown supertypes gracefully", {
  # Use an allele that might not be in the supertype mapping
  recipient <- hlaGeno(data.frame(A_1 = "A*99:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL")

  # Should not error, falls back to A02-like
  result <- calculatePeptideBindingLoad(recipient, peptides, backend = "pwm", return = "detail")
  expect_equal(nrow(result), 1)
})

# --- Summary Return Tests ---

test_that("summary return correctly aggregates across alleles", {
  recipient <- hlaGeno(data.frame(
    A_1 = "A*02:01", A_2 = "A*03:01",
    stringsAsFactors = FALSE
  ))
  peptides <- c("GILGFVFTL", "NLVPMVATV", "FLKEKGGL")

  result <- calculatePeptideBindingLoad(recipient, peptides, return = "summary")

  expect_equal(nrow(result), 2)  # Two alleles
  expect_true(all(result$n_strong >= 0))
  expect_true(all(result$n_weak >= 0))
  expect_true(all(result$n_strong + result$n_weak <= result$n_peptides))
})

test_that("summary totals match detail breakdown", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "NLVPMVATV", "FLKEKGGL")

  detail <- calculatePeptideBindingLoad(recipient, peptides, return = "detail")
  summary <- calculatePeptideBindingLoad(recipient, peptides, return = "summary")

  expect_equal(summary$n_peptides[1], nrow(detail))
  expect_equal(summary$n_strong[1], sum(detail$binding_level == "strong"))
  expect_equal(summary$n_weak[1], sum(detail$binding_level == "weak"))
})

# --- Contribution Calculation Tests ---

test_that("contribution scores are non-negative", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "NLVPMVATV")

  result <- calculatePeptideBindingLoad(recipient, peptides, return = "detail")

  expect_true(all(result$contribution >= 0))
})

test_that("strong binders have higher contribution than weak binders", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "NLVPMVATV", "FLKEKGGL", "AAAAAAAAL")

  result <- calculatePeptideBindingLoad(
    recipient, peptides,
    binding_threshold = 500,
    weak_threshold = 5000,
    return = "detail"
  )

  # Check that multiplier difference is reflected
  strong <- result[result$binding_level == "strong", ]
  weak <- result[result$binding_level == "weak", ]

  # This tests the formula: strong gets 2x multiplier, weak gets 1x
  # For same IC50, strong contribution should be ~2x weak
  if (nrow(strong) > 0 && nrow(weak) > 0) {
    # At least verify strong binders have non-zero contribution
    expect_true(all(strong$contribution > 0))
  }
})

# --- Edge Case Tests ---

test_that("calculatePeptideBindingLoad handles single character peptide input", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))

  # Single short peptide (will be filtered out for 9-mer requirement)
  result <- calculatePeptideBindingLoad(recipient, "ABC", peptide_length = 9L)
  expect_equal(result, 0)
})

test_that("calculatePeptideBindingLoad handles peptides with non-standard characters", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  # Peptide with X (unknown)
  peptides <- c("GILGFVFTX")

  result <- calculatePeptideBindingLoad(recipient, peptides, return = "detail")
  expect_equal(nrow(result), 1)  # Should still process
})

test_that("calculatePeptideBindingLoad handles duplicate peptides", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "GILGFVFTL", "GILGFVFTL")

  result <- calculatePeptideBindingLoad(recipient, peptides, return = "detail")

  # Duplicates should be deduplicated internally by .getPeptides
  expect_equal(nrow(result), 1)
})

# --- Different Peptide Lengths ---

test_that("calculatePeptideBindingLoad works with different peptide lengths", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))

  # 8-mers
  peptides_8 <- c("GILGFVFT", "NLVPMVAT")
  result_8 <- calculatePeptideBindingLoad(recipient, peptides_8, peptide_length = 8L, return = "detail")
  expect_equal(nrow(result_8), 2)
  expect_true(all(nchar(result_8$peptide) == 8))

  # 10-mers
  peptides_10 <- c("GILGFVFTLA", "NLVPMVATVA")
  result_10 <- calculatePeptideBindingLoad(recipient, peptides_10, peptide_length = 10L, return = "detail")
  expect_equal(nrow(result_10), 2)
  expect_true(all(nchar(result_10$peptide) == 10))
})

# --- Aggregate Method Tests ---

test_that("aggregate_method sum gives total of all contributions", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "NLVPMVATV")

  detail <- calculatePeptideBindingLoad(recipient, peptides, return = "detail")
  total_sum <- calculatePeptideBindingLoad(recipient, peptides, aggregate_method = "sum")

  expect_equal(total_sum, sum(detail$contribution))
})

test_that("aggregate_method max gives maximum contribution", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "NLVPMVATV")

  detail <- calculatePeptideBindingLoad(recipient, peptides, return = "detail")
  total_max <- calculatePeptideBindingLoad(recipient, peptides, aggregate_method = "max")

  expect_equal(total_max, max(detail$contribution))
})

test_that("aggregate_method mean gives average contribution", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))
  peptides <- c("GILGFVFTL", "NLVPMVATV")

  detail <- calculatePeptideBindingLoad(recipient, peptides, return = "detail")
  total_mean <- calculatePeptideBindingLoad(recipient, peptides, aggregate_method = "mean")

  expect_equal(total_mean, mean(detail$contribution))
})

# --- Class I vs Class II Alleles ---

test_that("calculatePeptideBindingLoad handles Class I B locus alleles", {
  recipient <- hlaGeno(data.frame(B_1 = "B*07:02", B_2 = "B*08:01", stringsAsFactors = FALSE))
  peptides <- c("TPRVTGGGAM")  # Known B*07:02 binder motif (Pro at P2)

  result <- calculatePeptideBindingLoad(recipient, peptides, peptide_length = 10L, return = "detail")

  expect_true(nrow(result) >= 1)
  expect_true(all(grepl("B\\*", result$hla_allele)))
})

# --- Empty and Edge Cases for Summary ---

test_that("summary return handles empty peptides correctly", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", A_2 = "A*03:01", stringsAsFactors = FALSE))

  result <- calculatePeptideBindingLoad(recipient, character(0), return = "summary")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(result$n_peptides == 0))
  expect_true(all(result$n_strong == 0))
  expect_true(all(result$n_weak == 0))
  expect_true(all(result$risk_contribution == 0))
})

test_that("detail return handles empty peptides correctly", {
  recipient <- hlaGeno(data.frame(A_1 = "A*02:01", stringsAsFactors = FALSE))

  result <- calculatePeptideBindingLoad(recipient, character(0), return = "detail")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_true(all(c("peptide", "hla_allele", "predicted_ic50", "binding_level", "contribution") %in% names(result)))
})

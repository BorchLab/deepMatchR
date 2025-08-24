# tests/testthat/test-spiDeconvolute.R

library(testthat)
library(data.table)

# No need to load external data, we will create it here.

test_that("spiDeconvolute returns correct structure with valid inputs", {
  sab_test <- data.table(antigen = "A1", NormalValue = 2000, class = "I")
  pra_test <- data.table(BeadID = 1, NormalValue = 1500)
  panel_test <- data.table(BeadID = 1, antigen = "A1")

  result <- spiDeconvolute(SAB = sab_test, PRA = pra_test, panelLong = panel_test)

  expect_true(is.list(result))
  expect_named(result, c("concordant", "sabOnly", "praOnly", "supportByBead", "summary", "argsUsed"))
  expect_true(is.data.frame(result$concordant))
  expect_true(is.data.frame(result$sabOnly))
  expect_true(is.data.frame(result$praOnly))
  expect_true(is.data.frame(result$supportByBead))
  expect_true(is.list(result$summary))
  expect_true(is.list(result$argsUsed))
})

test_that("spiDeconvolute partitions results correctly based on default cutoffs", {
  sab_test <- data.table(
    antigen = c("A1", "A2", "B7", "B8", "DR1"),
    NormalValue = c(2000, 500, 3000, 400, 4000),
    class = c("I", "I", "I", "I", "II")
  )
  pra_test <- data.table(
    BeadID = c(1, 2, 3, 4),
    NormalValue = c(1500, 400, 2000, 3000)
  )
  panel_test <- data.table(
    BeadID = c(1, 1, 2, 3, 3, 4),
    antigen = c("A1", "B8", "A2", "B7", "B8", "DR1")
  )

  # SAB cutoffs: Class I = 1500, Class II = 2500
  # Reactive: A1 (2000), B7 (3000), DR1 (4000)
  # Not reactive: A2 (500), B8 (400)

  # PRA cutoff: 1000
  # Positive beads: 1, 3, 4
  # Negative beads: 2

  # Support (positiveBeads >= 1):
  # A1: supported by Bead 1. Support = 1. -> Concordant
  # A2: supported by Bead 2 (neg). Support = 0. -> No category (correctly)
  # B7: supported by Bead 3. Support = 1. -> Concordant
  # B8: supported by Beads 1, 3. Support = 2. -> PRA-only
  # DR1: supported by Bead 4. Support = 1. -> Concordant

  result <- spiDeconvolute(SAB = sab_test, PRA = pra_test, panelLong = panel_test)

  expect_equal(sort(result$concordant$antigen), c("A1", "B7", "DR1"))
  expect_equal(nrow(result$sabOnly), 0)
  expect_equal(result$praOnly$antigen, "B8")
})

test_that("spiDeconvolute handles custom arguments (minSupportBeads, cutoffs)", {
  sab_test <- data.table(antigen = c("A1", "B8"), NormalValue = c(2000, 2500), class = "I")
  pra_test <- data.table(BeadID = c(1, 2), NormalValue = c(800, 900))
  panel_test <- data.table(BeadID = c(1, 2), antigen = c("A1", "B8"))

  args_list <- list(
    sabMfiCutoff = 1800,
    praMfiCutoff = 700,
    minSupportBeads = 2
  )

  # SAB cutoff = 1800. Reactive: A1, B8
  # PRA cutoff = 700. Positive beads: 1, 2
  # Support:
  # A1: supported by Bead 1. Support = 1.
  # B8: supported by Bead 2. Support = 1.
  # minSupportBeads = 2.
  # A1 is SAB reactive, but support (1) < minSupportBeads (2) -> sabOnly
  # B8 is SAB reactive, but support (1) < minSupportBeads (2) -> sabOnly

  result <- spiDeconvolute(SAB = sab_test, PRA = pra_test, panelLong = panel_test, args = args_list)

  expect_equal(nrow(result$concordant), 0)
  expect_equal(sort(result$sabOnly$antigen), c("A1", "B8"))
  expect_equal(nrow(result$praOnly), 0)
})

test_that("spiDeconvolute works with robust cutoff estimation", {
  sab_test <- data.table(antigen = paste0("A", 1:10),
                         NormalValue = c(rep(100, 8), 5000, 6000),
                         class = "I")
  pra_test <- data.table(BeadID = 1, NormalValue = 2000)
  panel_test <- data.table(BeadID = 1, antigen = "A10")

  result <- spiDeconvolute(SAB = sab_test, PRA = pra_test, panelLong = panel_test,
                           args = list(sabMfiCutoff = "robust"))

  # Robust cutoff is high, only A9 (5000) and A10 (6000) should be reactive.
  # From these two, only A10 is in the panel and supported by a positive bead.
  # So, A10 is concordant, A9 is sabOnly.
  expect_equal(result$concordant$antigen, "A10")
  expect_equal(result$sabOnly$antigen, "A9")
  expect_equal(nrow(result$praOnly), 0)
  # Check if the calculated cutoff is reasonable (should be > 1000)
  expect_true(result$argsUsed$sabMfiCutoff == "robust")
  expect_true(all(result$concordant$sab_cutoff > 1000))
})

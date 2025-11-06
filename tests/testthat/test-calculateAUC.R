#tests/testthat/test-calculateAUC.R
sab <- deepMatchR_example[[1]]

test_that("calculateAUC (eplet): returns ggplot and data, applies filters, computes norm_AUC", {
  # Plot path
  p <- calculateAUC(
    result_file   = sab,
    analysis_type = "eplet",
    cut_min = 250, cut_max = 5000, cut_step = 250,
    plot_results  = TRUE,
    label = FALSE
  )
  expect_true(ggplot2::is_ggplot(p))
  
  # Data path
  auc_dt <- calculateAUC(
    result_file   = sab,
    analysis_type = "eplet",
    cut_min = 250, cut_max = 5000, cut_step = 250,
    plot_results  = FALSE,
    label = FALSE
  )
  expect_s3_class(auc_dt, "data.frame")
  expect_true(all(c("eplet","AUC","norm_AUC","total_count","loci") %in% colnames(auc_dt)))
  expect_true(all(auc_dt$AUC >= 0))
  expect_true(all(auc_dt$norm_AUC >= 0 & auc_dt$norm_AUC <= 1))
})

test_that("calculateAUC (creg/serology) run and wrapper aliases dispatch", {
  
  
  # CREG
  creg_plot <- calculateAUC(sab, analysis_type = "creg", plot_results = TRUE, label = FALSE)
  expect_true(ggplot2::is.ggplot(creg_plot))
  creg_df   <- calculateAUC(sab, analysis_type = "creg", plot_results = FALSE, label = FALSE)
  expect_true(is.data.frame(creg_df))
  expect_true(all(c("CREG","AUC","norm_AUC","total_count","loci") %in% names(creg_df)))
  
  # Serology
  sero_plot <- calculateAUC(sab, analysis_type = "serology", plot_results = TRUE, label = FALSE)
  expect_true(ggplot2::is.ggplot(sero_plot))
  sero_df   <- calculateAUC(sab, analysis_type = "serology", plot_results = FALSE, label = FALSE)
  expect_true(is.data.frame(sero_df))
  expect_true(all(c("serology","AUC","norm_AUC","total_count","loci") %in% names(sero_df)))
  
  # Wrapper functions
  p1 <- epletAUC(sab, top_eplets = 5)
  expect_true(ggplot2::is_ggplot(p1))
  p2 <- cregAUC(sab)
  expect_true(ggplot2::is_ggplot(p2))
  p3 <- serologyAUC(sab)
  expect_true(ggplot2::is_ggplot(p3))
})

test_that("calculateAUC respects evidence_level filters and errors when filter empties set", {

  # Valid evidence subset
  expect_silent(
    calculateAUC(sab, analysis_type = "eplet",
                 evidence_level = c("A1","A2"),
                 plot_results = FALSE, label = FALSE)
  )
  
  # Evidence subset that will intentionally empty (no "Z9" in fixture)
  expect_error(
    calculateAUC(sab, analysis_type = "eplet",
                 evidence_level = "Z9",
                 plot_results = FALSE, label = FALSE),
    "did not produce any results"
  )
})

test_that("calculateAUC input validation", {
  
  expect_error(calculateAUC(sab, analysis_type = "wrong"), "must be one of")
  # Ensure top_eplets cropping doesn't error when fewer features exist
  expect_silent(
    calculateAUC(sab, analysis_type = "eplet", top_eplets = 100, label = FALSE)
  )
})

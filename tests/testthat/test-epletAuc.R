# test script for epletAUC.R - testcases are NOT comprehensive!

sab_data_example <- deepMatchR_example[[1]]

test_that("epletAUC() works with data frame input directly", {

  result_plot <- epletAUC(
    result_file = sab_data_example,
    plot_results = TRUE
  )

  expect_s3_class(result_plot, "ggplot")
})


test_that("epletAUC() returns data.frame with correct columns when plot_results=FALSE", {
  # Return a data.frame summarizing AUC
  result_df <- epletAUC(
    result_file = sab_data_example,
    plot_results = FALSE
  )

  expect_s3_class(result_df, "data.frame")
  expect_true(all(c("eplet", "AUC", "total_count") %in% colnames(result_df)))
})


test_that("epletAUC() throws error if missing SAB columns", {
  # Create a data frame missing a required column
  sab_data_incomplete <- data.table::as.data.table(sab_data_example)
  sab_data_incomplete[, BeadID := NULL]

  expect_error(
    epletAUC(
      result_file = sab_data_incomplete,
      plot_results = FALSE
    ),
    "Please ensure the SAB file includes the 'BeadID', 'SpecAbbr', 'Specificity', 'NormalValue'"
  )
})


test_that("epletAUC() processes numeric cut_min, cut_max, cut_step properly", {
  result_df <- epletAUC(
    result_file = sab_data_example,
    plot_results = FALSE,
    cut_min = 100,
    cut_max = 2000,
    cut_step = 500
  )

  expect_s3_class(result_df, "data.frame")
  # Must contain the usual columns
  expect_true(all(c("eplet", "AUC", "total_count") %in% names(result_df)))
})

test_that("epletAUC() can handle different evidence_level inputs", {
  expect_error(epletAUC(
    result_file     = sab_data_example,
    plot_results    = FALSE,
    evidence_level  = "Nonexistent_Level"
  ))
})


test_that("epletAUC() final tibble has expected numeric values for AUC", {
  result_df <- epletAUC(
    result_file  = sab_data_example,
    plot_results = FALSE
  )

  expect_type(result_df$AUC, "double")
  expect_type(result_df$total_count, "integer")
})

# --- Tests for cregAUC ---

test_that("cregAUC() returns ggplot when plot_results=TRUE", {
  result_plot <- cregAUC(
    result_file = sab_data_example,
    plot_results = TRUE
  )

  expect_s3_class(result_plot, "ggplot")
})

test_that("cregAUC() returns data.frame when plot_results=FALSE", {
  result_df <- cregAUC(
    result_file = sab_data_example,
    plot_results = FALSE
  )

  expect_s3_class(result_df, "data.frame")
  expect_true("CREG" %in% names(result_df))
  expect_true("AUC" %in% names(result_df))
})

test_that("cregAUC() respects creg_filter parameter", {
  result_strict <- cregAUC(
    result_file = sab_data_example,
    creg_filter = 10,
    plot_results = FALSE
  )

  result_loose <- cregAUC(
    result_file = sab_data_example,
    creg_filter = 1,
    plot_results = FALSE
  )

  # Stricter filter should return fewer or equal rows
  expect_true(nrow(result_strict) <= nrow(result_loose))
})

# --- Tests for serologyAUC ---

test_that("serologyAUC() returns ggplot when plot_results=TRUE", {
  result_plot <- serologyAUC(
    result_file = sab_data_example,
    plot_results = TRUE
  )

  expect_s3_class(result_plot, "ggplot")
})

test_that("serologyAUC() returns data.frame when plot_results=FALSE", {
  result_df <- serologyAUC(
    result_file = sab_data_example,
    plot_results = FALSE
  )

  expect_s3_class(result_df, "data.frame")
  expect_true("serology" %in% names(result_df))
  expect_true("AUC" %in% names(result_df))
})

test_that("serologyAUC() respects serology_filter parameter", {
  result_strict <- serologyAUC(
    result_file = sab_data_example,
    serology_filter = 10,
    plot_results = FALSE
  )

  result_loose <- serologyAUC(
    result_file = sab_data_example,
    serology_filter = 1,
    plot_results = FALSE
  )

  # Stricter filter should return fewer or equal rows
  expect_true(nrow(result_strict) <= nrow(result_loose))
})

# --- Tests for calculateAUC base function ---

test_that("calculateAUC() errors on invalid analysis_type", {
  expect_error(
    calculateAUC(sab_data_example, analysis_type = "invalid"),
    "analysis_type"
  )
})

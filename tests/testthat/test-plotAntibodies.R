# tests/testthat/test-plotAntibodies.R 

test_that("plotAntibodies returns a ggplot object for SAB data (add_table = FALSE)", {
  sab_data <- deepMatchR_example[[1]]

  p <- plotAntibodies(result_file = sab_data, type = "SAB", add_table = FALSE)
  expect_s3_class(p, "ggplot")
  expect_false("patchwork" %in% class(p))
})

test_that("plotAntibodies returns a patchwork object for SAB data when add_table = TRUE", {
  sab_data <- deepMatchR_example[[1]]

  p <- plotAntibodies(result_file = sab_data, type = "SAB", add_table = TRUE)
  expect_true("patchwork" %in% class(p))
})

test_that("plotAntibodies correctly processes file path input for SAB data", {
  sab_data <- deepMatchR_example[[1]]

  temp_file <- tempfile(fileext = ".csv")
  write.csv(sab_data, temp_file, row.names = FALSE)

  p <- plotAntibodies(result_file = temp_file, type = "SAB", add_table = FALSE)
  expect_s3_class(p, "ggplot")

  unlink(temp_file)
})

test_that("plotAntibodies highlights antigens for SAB data", {
  sab_data <- deepMatchR_example[[1]]

  expect_error(
    plotAntibodies(result_file = sab_data, type = "SAB", highlight_antigen = "NonExistentAntigen"),
    regexp = "highlight_antigen selection is not within the resulting data.frame"
  )

  valid_antigen <- "A1"
  p <- plotAntibodies(result_file = sab_data, type = "SAB", highlight_antigen = valid_antigen)
  expect_s3_class(p, "ggplot")
})

test_that("plotAntibodies returns a ggplot object for PRA data (add_table = FALSE)", {
  pra_data <- deepMatchR_example[[3]]

  p <- plotAntibodies(result_file = pra_data, type = "PRA", class = "I", add_table = FALSE)
  expect_s3_class(p, "ggplot")
  expect_false("patchwork" %in% class(p))
})

test_that("plotAntibodies returns a patchwork object for PRA data when add_table = TRUE", {
  pra_data <- deepMatchR_example[[3]]

  p <- plotAntibodies(result_file = pra_data, type = "PRA", class = "I", add_table = TRUE)
  expect_true("patchwork" %in% class(p))
})

test_that("plotAntibodies returns a ggplot object for trend plot", {
  sab_data_list <- list(
    "01/01/2023" = deepMatchR_example[[1]],
    "02/01/2023" = deepMatchR_example[[1]]
  )

  p <- plotAntibodies(result_file = sab_data_list, type = "SAB", plot_trend = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("plotAntibodies trend plot errors with incorrect input", {
  sab_data <- deepMatchR_example[[1]]

  expect_error(
    plotAntibodies(result_file = sab_data, type = "SAB", plot_trend = TRUE),
    regexp = "For trend plots, 'result_file' must be a named list of data frames."
  )
})

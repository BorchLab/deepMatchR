# test script for plotEplets.R - testcases are NOT comprehensive!

test_that("plotEplets returns a ggplot object for treemap", {
  p <- plotEplets(result_file = deepMatchR_example[[1]],
                  plot.type = "treemap",
                  cutoff = 2000,
                  evidence_level = c("A1", "A2"),
                  percPos_filter = 0.4,
                  palette = "inferno")
  expect_true(inherits(p, "ggplot"))
})

test_that("plotEplets returns a ggplot object for bar", {
  p <- plotEplets(result_file = deepMatchR_example[[1]],
                  plot.type = "bar",
                  cutoff = 2000,
                  evidence_level = c("A1", "A2"),
                  percPos_filter = 0.4,
                  top.eplets = 20,
                  palette = "inferno")
  expect_true(inherits(p, "ggplot"))
})

test_that("plotEplets returns a ggplot object for AUC", {
  p <- plotEplets(result_file = deepMatchR_example[[1]],
                  plot.type = "AUC",
                  percPos_filter = 0.4,
                  cut_min = 250,
                  cut_max = 10000,
                  cut_step = 250,
                  top.eplets = 20,
                  palette = "inferno")
  expect_true(inherits(p, "ggplot"))
})

test_that("plotEplets errors with an invalid plot type", {
  expect_error(
    plotEplets(result_file = deepMatchR_example[[1]], plot.type = "invalid"),
    regexp = "should be one of"
  )
})

test_that("plotEplets handles file path input", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(deepMatchR_example[[1]], tmp, row.names = FALSE)
  
  p <- plotEplets(result_file = tmp,
                  plot.type = "bar",
                  cutoff = 2000,
                  evidence_level = c("A1", "A2"),
                  percPos_filter = 0.4,
                  top.eplets = 20,
                  palette = "inferno")
  expect_true(inherits(p, "ggplot"))
  unlink(tmp)
})


test_that("plotEplets works with evidence_level NULL", {
  expect_silent(
    plotEplets(result_file = deepMatchR_example[[1]],
               evidence_level = NULL,
               plot.type = "bar",
               cutoff = 2000,
               percPos_filter = 0.4,
               top.eplets = 20,
               palette = "inferno")
  )
})


test_that("plotEplets errors if required SAB columns are missing", {
  # Create a copy of the example data with a required column removed.
  bad_data <- deepMatchR_example[[1]]
  bad_data$NormalValue <- NULL
  
  expect_error(
    plotEplets(result_file = bad_data,
               plot.type = "treemap",
               cutoff = 2000,
               evidence_level = c("A1", "A2"),
               percPos_filter = 0.4,
               palette = "inferno"),
    regexp = "Please ensure the SAB file includes the"
  )
})

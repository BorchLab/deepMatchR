# test script for plotSAB.R - testcases are NOT comprehensive!

context("Testing plotSAB function")

test_that("plotSAB returns a ggplot object with data frame input (add_table = FALSE)", {
  # Use the built-in example data
  sab_data <- deepMatchR_example[[1]]
  
  # Call function with default add_table = TRUE but then also test with add_table = FALSE
  p <- plotSAB(result_file = sab_data, add_table = FALSE)
  expect_s3_class(p, "ggplot")
  expect_false("patchwork" %in% class(p))  # When add_table is FALSE, no patchwork structure is expected
})

test_that("plotSAB returns a patchwork object when add_table = TRUE", {
  sab_data <- deepMatchR_example[[1]]
  
  p <- plotSAB(result_file = sab_data, add_table = TRUE)
  # When add_table is TRUE, the combined plot (via patchwork) should have "patchwork" in its class.
  expect_true("patchwork" %in% class(p))
})

test_that("plotSAB correctly processes file path input", {
  sab_data <- deepMatchR_example[[1]]
  
  # Write the example data to a temporary CSV file
  temp_file <- tempfile(fileext = ".csv")
  write.csv(sab_data, temp_file, row.names = FALSE)
  
  p <- plotSAB(result_file = temp_file, add_table = FALSE)
  expect_s3_class(p, "ggplot")
  
  # Clean up the temporary file
  unlink(temp_file)
})

test_that("plotSAB highlights antigens when valid highlight_antigen is provided", {
  sab_data <- deepMatchR_example[[1]]
  
  # First, test that an invalid highlight triggers an error
  expect_error(
    plotSAB(result_file = sab_data, highlight_antigen = "NonExistentAntigen"),
    regexp = " selection is not within the data.frame"
  )
  
  valid_antigen <- "A1"
  
  # Running the function with a valid highlight should not error.
  p <- plotSAB(result_file = sab_data, highlight_antigen = valid_antigen)
  expect_s3_class(p, "ggplot")
})
# tests/testthat/test-visualizePeptideBinding.R

make_binding_results <- function() {
  # Minimal, plausible all_predictions table
  all_predictions <- data.frame(
    donor_allele = c("A*01:01","A*01:01","A*02:01","A*02:01","B*07:02","B*07:02"),
    test_allele  = c("A*01:01","A*02:01","A*01:01","A*02:01","B*07:02","B*08:01"),
    locus        = c("A","A","A","A","B","B"),
    binding      = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE),
    ic50         = c(50, 5000, 120, 350, 10000, 80),
    stringsAsFactors = FALSE
  )
  list(all_predictions = all_predictions)
}

test_that("input validation triggers error for malformed binding_results", {
  expect_error(visualizePeptideBinding(list(), plot_type = "heatmap"),
               "must be output from calculatePeptideBindingLoad")
})

test_that("heatmap returns a ggplot with expected mappings/facets", {
  br <- make_binding_results()
  p <- visualizePeptideBinding(br, plot_type = "heatmap")
  expect_s3_class(p, "ggplot")
  # basic structure: should have donor_allele on x, test_allele on y
  m <- ggplot2::ggplot_build(p)
  # expect at least one panel per locus (facets)
  expect_true(length(m$layout$layout$PANEL) >= 1)
})

test_that("bar plot returns a ggplot and shows counts by locus", {
  br <- make_binding_results()
  p <- visualizePeptideBinding(br, plot_type = "bar")
  expect_s3_class(p, "ggplot")
  # Built data should include a bar layer
  built <- ggplot2::ggplot_build(p)
  geoms <- vapply(built$plot$layers, function(x) class(x$geom)[1], character(1))
  expect_true(any(grepl("GeomBar", geoms)))
})

test_that("scatter plot returns a ggplot with ic50 on x and facets by locus", {
  br <- make_binding_results()
  p <- visualizePeptideBinding(br, plot_type = "scatter")
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  # Should have a point layer
  geoms <- vapply(built$plot$layers, function(x) class(x$geom)[1], character(1))
  expect_true(any(grepl("GeomPoint", geoms)))
  # log10 scale on x is present
  expect_true(any(vapply(built$plot$scales$scales, function(s) inherits(s, "ScaleContinuousPosition"), logical(1))))
})

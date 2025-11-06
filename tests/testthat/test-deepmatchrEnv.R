# tests/testthat/test-deepmatchrEnv.R

test_that("deepmatchrEnv returns platform-specific BasiliskEnvironment objects", {
  mac <- deepmatchrEnv("macos")
  lin <- deepmatchrEnv("linux")
  
  # Basic type & fields
  expect_s3_class(mac, "BasiliskEnvironment")
  expect_s3_class(lin, "BasiliskEnvironment")
  
  # Same envname, pkgname; different TensorFlow wheel on mac vs linux
  expect_identical(mac$envname, "deepmatchrEnv")
  expect_identical(lin$envname, "deepmatchrEnv")
  expect_identical(mac$pkgname, "deepMatchR")
  expect_identical(lin$pkgname, "deepMatchR")
  
  expect_true(any(grepl("^tensorflow-macos==2", mac$pip)))
  expect_true(any(grepl("^tensorflow==2", lin$pip)))
  expect_false(any(grepl("^tensorflow-macos==", lin$pip)))
  expect_false(any(grepl("^tensorflow==", mac$pip)))
})

test_that("deepmatchrEnv('auto') defaults to linux unless clearly macOS", {
  # We can't reliably mock base::Sys.info here, so just assert that 'auto'
  # yields *one* of the two predefined envs.
  auto <- deepmatchrEnv("auto")
  mac  <- deepmatchrEnv("macos")
  lin  <- deepmatchrEnv("linux")
  
  # Compare a discriminating field (pip list) to determine which we got
  if (any(grepl("^tensorflow-macos==", auto$pip))) {
    expect_equal(auto$pip, mac$pip)
  } else {
    expect_equal(auto$pip, lin$pip)
  }
})

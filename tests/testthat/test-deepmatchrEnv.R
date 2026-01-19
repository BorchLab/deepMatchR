# tests/testthat/test-deepmatchrEnv.R

# Adjust if your package name differs:
.pkg <- "deepMatchR"

# Helpers to fetch the internal S4 env objects from your package namespace
get_linux_env  <- function() getFromNamespace(".deepmatchrEnv_linux", .pkg)
get_macos_env  <- function() getFromNamespace(".deepmatchrEnv_macos", .pkg)

sysname <- tryCatch({
  info <- Sys.info()
  if (!is.null(info) && length(info) && !is.null(info[["sysname"]])) {
    tolower(info[["sysname"]])
  } else {
    tolower(.Platform$OS.type)
  }
}, error = function(e) {
  tolower(.Platform$OS.type)
})

test_that("explicit platform selection returns expected BasiliskEnvironment S4 objects", {
  out_lin <- deepmatchrEnv("linux")
  out_mac <- deepmatchrEnv("macos")
  
  # S4 class checks
  expect_s4_class(out_lin, "BasiliskEnvironment")
  expect_s4_class(out_mac, "BasiliskEnvironment")
  
  # envname slot is the one you defined in your code ("deepmatchrEnv")
  expect_identical(out_lin@envname, "deepmatchrEnv_v2")
  expect_identical(out_mac@envname, "deepmatchrEnv_v2")
})

test_that("match.arg validates unsupported platform", {
  expect_error(deepmatchrEnv("windows"), "arg")
})

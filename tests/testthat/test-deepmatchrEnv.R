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

# --- Auto Platform Detection Tests ---

test_that("auto platform detection returns valid BasiliskEnvironment", {
  out_auto <- deepmatchrEnv("auto")
  expect_s4_class(out_auto, "BasiliskEnvironment")
  expect_identical(out_auto@envname, "deepmatchrEnv_v2")
})

test_that("auto platform matches explicit selection for current OS", {
  out_auto <- deepmatchrEnv("auto")

  if (grepl("darwin|mac", sysname)) {
    out_explicit <- deepmatchrEnv("macos")
  } else {
    out_explicit <- deepmatchrEnv("linux")
  }

  expect_identical(out_auto@envname, out_explicit@envname)
  expect_identical(out_auto@pkgname, out_explicit@pkgname)
})

# --- Environment Configuration Tests ---

test_that("linux environment has correct package configuration", {
  out_lin <- deepmatchrEnv("linux")

  expect_identical(out_lin@pkgname, "deepMatchR")
  # Check pip packages are specified (the slots vary by basilisk version)
  # At minimum, envname and pkgname should be correct
})

test_that("macos environment has correct package configuration", {
  out_mac <- deepmatchrEnv("macos")

  expect_identical(out_mac@pkgname, "deepMatchR")
})

# --- Default Argument Tests ---

test_that("deepmatchrEnv with no arguments uses auto", {
  out_default <- deepmatchrEnv()
  out_auto <- deepmatchrEnv("auto")

  expect_identical(out_default@envname, out_auto@envname)
})

# --- Internal Environment Object Tests ---

test_that("internal linux environment object exists and is valid", {
  env_lin <- get_linux_env()
  expect_s4_class(env_lin, "BasiliskEnvironment")
  expect_identical(env_lin@envname, "deepmatchrEnv_v2")
  expect_identical(env_lin@pkgname, "deepMatchR")
})

test_that("internal macos environment object exists and is valid", {
  env_mac <- get_macos_env()
  expect_s4_class(env_mac, "BasiliskEnvironment")
  expect_identical(env_mac@envname, "deepmatchrEnv_v2")
  expect_identical(env_mac@pkgname, "deepMatchR")
})

# --- Idempotency Tests ---

test_that("repeated calls return identical environments", {
  out1 <- deepmatchrEnv("linux")
  out2 <- deepmatchrEnv("linux")

  expect_identical(out1@envname, out2@envname)
  expect_identical(out1@pkgname, out2@pkgname)
})

test_that("environment names are consistent across platforms", {
  out_lin <- deepmatchrEnv("linux")
  out_mac <- deepmatchrEnv("macos")

  # Both should have the same environment name for consistency
  expect_identical(out_lin@envname, out_mac@envname)
})

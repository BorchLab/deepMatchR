# tests/testthat/test-clearSequenceCache.R
library(withr)

test_that("clearing non-existent filesystem cache emits message and does nothing", {
  td <- file.path(tempdir(), paste0("no_such_cache_", as.integer(runif(1, 1e6, 2e6))))
  expect_false(dir.exists(td))
  expect_message(clearSequenceCache(cache_dir = td), "does not exist")
  expect_false(dir.exists(td))
})

test_that("clearing whole filesystem cache removes directory", {
  td <- local_tempdir()
  # Create some dummy files
  f1 <- file.path(td, "foo_PROT_bar.rds")
  f2 <- file.path(td, "baz_NUC_qux.rds")
  f3 <- file.path(td, "other_anything.txt")
  writeLines("x", f1)
  writeLines("x", f2)
  writeLines("x", f3)
  expect_true(file.exists(f1) && file.exists(f2) && file.exists(f3))
  
  expect_message(clearSequenceCache(cache_dir = td), "Cleared filesystem cache")
  # Directory should be gone
  expect_false(dir.exists(td))
})

test_that("type-specific clearing removes only matching files", {
  td <- local_tempdir()
  # PROT files
  fP1 <- file.path(td, "a_PROT_1.rds")
  fP2 <- file.path(td, "b_PROT_2.rds")
  # NUC files
  fN1 <- file.path(td, "a_NUC_1.rds")
  fN2 <- file.path(td, "b_NUC_2.rds")
  # other files
  fO  <- file.path(td, "unrelated.rds")
  
  for (f in c(fP1,fP2,fN1,fN2,fO)) writeLines("x", f)
  expect_true(all(file.exists(c(fP1,fP2,fN1,fN2,fO))))
  
  # Clear only PROT
  expect_message(clearSequenceCache(cache_dir = td, type = "PROT"), "Cleared 2 PROT cache files")
  expect_false(file.exists(fP1))
  expect_false(file.exists(fP2))
  expect_true(file.exists(fN1))
  expect_true(file.exists(fN2))
  expect_true(file.exists(fO))
  
  # Clear only NUC
  expect_message(clearSequenceCache(cache_dir = td, type = "NUC"), "Cleared 2 NUC cache files")
  expect_false(file.exists(fN1))
  expect_false(file.exists(fN2))
  expect_true(file.exists(fO))
  
  # Final clear-all
  expect_message(clearSequenceCache(cache_dir = td), "Cleared filesystem cache")
  expect_false(dir.exists(td))
})

test_that("in-memory clear path emits informative message and returns invisibly", {
  expect_message(clearSequenceCache(), "In-memory cache clearing", fixed = TRUE)
  expect_invisible(clearSequenceCache())
})

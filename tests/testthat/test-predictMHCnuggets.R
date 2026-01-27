# tests/testthat/test-predictMHCnuggets.R

testthat::local_edition(3)
testthat::skip_if_not_installed("basilisk")
testthat::skip_if_not_installed("reticulate")

# Test Windows error handling first (before skipping on Windows)
test_that("predictMHCnuggets errors on Windows", {
 testthat::skip_if_not(.Platform$OS.type == "windows", "Only runs on Windows")
 expect_error(
    predictMHCnuggets(peptides = "SIINFEKL", allele = "A*02:01"),
    "not supported on Windows"
  )
})

# Skip remaining tests on Windows since MHCnuggets is not supported
testthat::skip_on_os("windows")

# Capture arguments passed through basiliskRun for assertions
.observed <- new.env(parent = emptyenv())

# Mocks
mock_basiliskStart <- function(env) "proc-token"
mock_basiliskStop  <- function(proc) invisible(TRUE)

# By default, basiliskRun will write a CSV to out_path and return TRUE.
mock_basiliskRun_success <- function(proc, fun, ..., out_path) {
  # Capture forwarded args
  dots <- list(...)
  .observed$mhc      <- dots$mhc
  .observed$class    <- dots$cls
  .observed$model    <- dots$model
  .observed$out_path <- out_path
  
  # Simulate Python produced CSV
  df <- data.frame(
    peptide = c("SIINFEKL","LLFGYPVYV"),
    ic50    = c(23.5, 1234.0),
    rank    = c(0.01, 0.25),
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, out_path, row.names = FALSE)
  TRUE
}

mock_basiliskRun_affinity <- function(proc, fun, ..., out_path) {
  # Write CSV with 'affinity' instead of 'ic50'
  df <- data.frame(
    peptide  = c("AAA","BBB"),
    affinity = c(500, 20000),
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, out_path, row.names = FALSE)
  TRUE
}

mock_basiliskRun_fail <- function(proc, fun, ..., out_path) {
  # Do not write any file, and return FALSE
  FALSE
}

# Reticulate mocks (no-ops)
mock_py_run_string <- function(...) invisible(NULL)
mock_import        <- function(...) structure(list(predict = function(...) NULL), class = "pyobj")

# Helper macro: run a body with both package namespaces mocked
with_basilisk_and_reticulate <- function(body) {
  testthat::with_mocked_bindings(
    .package = "basilisk",
    basiliskStart = mock_basiliskStart,
    basiliskRun   = mock_basiliskRun_success,  # you can swap per-test
    basiliskStop  = mock_basiliskStop,
    {
      testthat::with_mocked_bindings(
        .package = "reticulate",
        py_run_string = mock_py_run_string,
        import        = mock_import,
        { force(body) }
      )
    }
  )
}

test_that("predictMHCnuggets returns a data.frame with expected columns", {
  testthat::with_mocked_bindings(
    .package = "basilisk",
    basiliskStart = mock_basiliskStart,
    basiliskRun   = mock_basiliskRun_success,
    basiliskStop  = mock_basiliskStop,
    {
      testthat::with_mocked_bindings(
        .package = "reticulate",
        py_run_string = mock_py_run_string,
        import        = mock_import,
        {
          res <- predictMHCnuggets(
            peptides   = c("SIINFEKL","LLFGYPVYV"),
            allele     = "A*02:01",
            mhc_class  = "I",
            hla_env    = deepmatchrEnv("linux"),
            rank_output = TRUE
          )
          expect_true(is.data.frame(res))
          expect_true(all(c("peptide","ic50") %in% names(res)))
          expect_true(is.numeric(res$ic50))
          expect_true(all(res$peptide %in% c("SIINFEKL","LLFGYPVYV")))
        }
      )
    }
  )
})

test_that("allele normalization feeds normalized mhc to basiliskRun (class I and II)", {
  testthat::with_mocked_bindings(
    .package = "basilisk",
    basiliskStart = mock_basiliskStart,
    basiliskRun   = mock_basiliskRun_success,
    basiliskStop  = mock_basiliskStop,
    {
      testthat::with_mocked_bindings(
        .package = "reticulate",
        py_run_string = mock_py_run_string,
        import        = mock_import,
        {
          # Class I: "A0201" -> "HLA-A02:01"
          invisible(predictMHCnuggets(
            peptides  = "SIINFEKL",
            allele    = "A0201",
            mhc_class = "I",
            hla_env   = deepmatchrEnv("linux")
          ))
          expect_identical(.observed$mhc, "HLA-A02:01")
          expect_identical(.observed$class, "I")
          
          # Class II: "DRB1*01:01" -> "HLA-DRB101:01"
          invisible(predictMHCnuggets(
            peptides  = "AAAABBBBCCCC",
            allele    = "DRB1*01:01",
            mhc_class = "II",
            hla_env   = deepmatchrEnv("linux")
          ))
          expect_identical(.observed$mhc, "HLA-DRB101:01")
          expect_identical(.observed$class, "II")
        }
      )
    }
  )
})

test_that("predictMHCnuggets coerces 'affinity' to 'ic50' when needed", {
  testthat::with_mocked_bindings(
    .package = "basilisk",
    basiliskStart = mock_basiliskStart,
    basiliskRun   = mock_basiliskRun_affinity,
    basiliskStop  = mock_basiliskStop,
    {
      testthat::with_mocked_bindings(
        .package = "reticulate",
        py_run_string = mock_py_run_string,
        import        = mock_import,
        {
          res <- predictMHCnuggets(
            peptides  = c("AAA","BBB"),
            allele    = "A*02:01",
            mhc_class = "I",
            hla_env   = deepmatchrEnv("linux")
          )
          expect_true("ic50" %in% names(res))
          expect_type(res$ic50, "double")
        }
      )
    }
  )
})

test_that("predictMHCnuggets errors when no CSV is produced", {
  testthat::with_mocked_bindings(
    .package = "basilisk",
    basiliskStart = mock_basiliskStart,
    basiliskRun   = mock_basiliskRun_fail,
    basiliskStop  = mock_basiliskStop,
    {
      testthat::with_mocked_bindings(
        .package = "reticulate",
        py_run_string = mock_py_run_string,
        import        = mock_import,
        {
          expect_error(
            predictMHCnuggets(peptides = "SIINFEKL", allele = "A*02:01", hla_env = deepmatchrEnv("linux")),
            "did not produce an output file"
          )
        }
      )
    }
  )
})

test_that("predictMHCnuggets validates inputs and returns empty df for length-0 peptides", {
  testthat::with_mocked_bindings(
    .package = "basilisk",
    basiliskStart = mock_basiliskStart,
    basiliskRun   = mock_basiliskRun_success,
    basiliskStop  = mock_basiliskStop,
    {
      testthat::with_mocked_bindings(
        .package = "reticulate",
        py_run_string = mock_py_run_string,
        import        = mock_import,
        {
          out <- predictMHCnuggets(character(), allele = "A*02:01", hla_env = deepmatchrEnv("linux"))
          expect_identical(nrow(out), 0L)
          expect_identical(names(out), "peptide")
          
          expect_error(predictMHCnuggets(42, "A*02:01"), "peptides")
          expect_error(predictMHCnuggets("AAA", c("A*02:01","A*01:01")), "single string")
          expect_error(predictMHCnuggets("AAA", "A*02:01", mhc_class = "III"), "must be 'I' or 'II'")
        }
      )
    }
  )
})

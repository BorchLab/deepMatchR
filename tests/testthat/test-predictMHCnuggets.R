# tests/testthat/test-predictMHCnuggets.R

test_that("predictMHCnuggets wrapper works", {
  
  # Mock basiliskRun
  mock_basiliskRun <- function(proc, fun, peptides_path, out_path, ...) {
    # 1. Check that inputs are correct
    expect_true(file.exists(peptides_path))
    peps <- readLines(peptides_path)
    expect_equal(peps, c("SIINFEKL", "LLFGYPVYV"))
    
    # 2. Check allele normalization
    # The `...` args are passed to `fun`, but `fun` is run by basilisk.
    # We can grab them from the parent environment
    run_args <- list(...)
    expect_equal(run_args$mhc, "HLA-A02:01")
    expect_equal(run_args$cls, "I")
    
    # 3. Create the dummy output file
    dummy_df <- data.frame(
      peptide = peps,
      ic50 = c(10.5, 200.0),
      rank = c(0.5, 2.0)
    )
    utils::write.csv(dummy_df, file = out_path, row.names = FALSE)
    
    # 4. Return TRUE (success)
    return(TRUE)
  }
  
  # Mock basiliskStart
  mock_basiliskStart <- function(...) { return(list(pid = 123)) } # Return a dummy process
  mock_basiliskStop <- function(...) { invisible(NULL) }
  
  testthat::with_mocked_bindings(
    basilisk::basiliskRun = mock_basiliskRun,
    basilisk::basiliskStart = mock_basiliskStart,
    basilisk::basiliskStop = mock_basiliskStop,
    .deepmatchrEnv_linux = list(), # Mock env
    deepmatchrEnv = function(...) list(), # Mock env
    {
      res <- predictMHCnuggets(
        peptides = c("SIINFEKL", "LLFGYPVYV"),
        allele = "A*02:01", # Test normalization
        mhc_class = "I",
        rank_output = TRUE,
        hla_env = list() # Pass dummy env
      )
      
      expect_s3_class(res, "data.frame")
      expect_equal(nrow(res), 2)
      expect_equal(res$ic50, c(10.5, 200.0))
      expect_equal(res$rank, c(0.5, 2.0))
    }
  )
})
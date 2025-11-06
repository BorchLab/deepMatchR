# tests/testthat/test-calculatePeptideBindingLoad.R

# Mock sequence DB
mock_seq_db <- c(
  "A*01:01" = "YFAMYGEKVAHTHVDTLYVRYHY",
  "A*02:01" = "YFDMYGEKVAHTHVDTLYVRFHY",
  "B*07:02" = "GSHSMRYFYTAMSRPGRGEPRFI",
  "B*08:01" = "GSHSMRYFYTAMSRPGRGEPRFI"
)

# Mocks
mock_validate <- function(x) TRUE
mock_batchGetSequences <- function(alleles, type = "PROT", ...) {
  miss <- setdiff(unique(alleles), names(mock_seq_db))
  if (length(miss)) stop("Unknown allele(s): ", paste(miss, collapse = ", "))
  mock_seq_db[unique(alleles)]
}
mock_predict <- function(peptides, allele, mhc_class = "I", ic50_threshold = 500, ...) {
  # Return a trivial binding table with rank + ic50 so downstream merges work
  data.frame(
    peptide = peptides,
    ic50 = rep(1000, length(peptides)),
    rank = rep(0.5, length(peptides)),
    stringsAsFactors = FALSE
  )
}

test_that("no shared loci -> error", {
  rgeno <- hlaGeno(data.frame(C_1 = "C*01:02", C_2 = "C*02:02", check.names = FALSE))
  dgeno <- hlaGeno(data.frame(A_1 = "A*01:01", A_2 = "A*02:01", check.names = FALSE))
  
  testthat::with_mocked_bindings(
    validateHlaGeno   = mock_validate,
    batchGetSequences = mock_batchGetSequences,
    predictMHCnuggets = mock_predict,
    {
      expect_error(calculatePeptideBindingLoad(rgeno, dgeno), "No shared loci")
    }
  )
})



test_that("END-TO-END: shapes for summary/per_locus/detailed", {

  rgeno <- hlaGeno(data.frame(A_1 = "A*01:01", A_2 = "A*02:01", B_1 = "B*07:02", B_2 = "B*08:01", check.names = FALSE))
  dgeno <- hlaGeno(data.frame(A_1 = "A*02:01", A_2 = "A*01:01", B_1 = "B*08:01", B_2 = "B*07:02", check.names = FALSE))
  
  testthat::with_mocked_bindings(
    validateHlaGeno   = mock_validate,
    batchGetSequences = mock_batchGetSequences,
    predictMHCnuggets = mock_predict,
    {
      s <- calculatePeptideBindingLoad(rgeno, dgeno, return = "summary", parallel = FALSE, mhc_class = "I")
      expect_true(is.data.frame(s))
      expect_true(all(c("total_mismatched_peptides","binding_peptides","binding_percentage","ic50_threshold","mhc_class") %in% names(s)))
      
      p <- calculatePeptideBindingLoad(rgeno, dgeno, return = "per_locus", parallel = FALSE)
      expect_true(is.data.frame(p))
      expect_true(all(c("locus","total_peptides","binding_peptides","binding_percentage","mean_ic50_binders") %in% names(p)))
      
      d <- calculatePeptideBindingLoad(rgeno, dgeno, return = "detailed", parallel = FALSE)
      expect_true(is.list(d))
      expect_true(is.data.frame(d$summary))
      expect_true(is.list(d$per_locus))
      expect_true(is.data.frame(d$all_predictions))
    }
  )
})

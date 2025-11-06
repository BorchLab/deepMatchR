#' Compare Mismatch Load vs Peptide Binding Load
#'
#' @description
#' Performs correlation analysis between traditional mismatch metrics and
#' peptide binding predictions
#'
#' @param recipient_geno Recipient genotype
#' @param donor_geno Donor genotype
#' @param loci Loci to analyze
#' @param mhc_class MHC class for peptide prediction
#'
#' @return Data frame with comparison metrics
#'
#' @export
compareLoadMetrics <- function(recipient_geno, 
                               donor_geno, 
                               loci = NULL, 
                               mhc_class = "I") {
  
  # Calculate mismatch load
  mismatch_load <- calculateMismatchLoadFast(
    recipient_geno, donor_geno, 
    loci = loci, 
    return = "per_locus"
  )
  
  # Calculate eplet load
  eplet_load <- calculateEpletLoad(
    recipient_geno, donor_geno,
    loci = loci,
    return = "per_locus"
  )
  
  # Calculate peptide binding load
  peptide_load <- calculatePeptideBindingLoad(
    recipient_geno, donor_geno,
    loci = loci,
    mhc_class = mhc_class,
    return = "per_locus"
  )
  
  # Merge results
  result <- merge(mismatch_load, eplet_load, by = "locus")
  result <- merge(result, peptide_load[, c("locus", "binding_peptides", "binding_percentage")], 
                  by = "locus")
  
  # Calculate correlations
  result$mismatch_eplet_ratio <- result$mismatch_load / (result$eplet_load + 1)
  result$binding_per_mismatch <- result$binding_peptides / (result$mismatch_load + 1)
  
  return(result)
}

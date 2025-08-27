#' Get Amino Acid Sequence for an HLA Allele
#'
#' @description
#' This function retrieves the amino acid sequence for a given HLA allele
#' name using the `immReferent` package.
#'
#' @param allele_name A character string representing the HLA allele name
#'   (e.g., "A*01:01").
#'
#' @return A character string representing the amino acid sequence.
#'
#' @importFrom immReferent getIMGT
#' @export
getAlleleSequence <- function(allele_name) {
  # Get all HLA protein sequences
  hla_sequences <- immReferent::getIMGT(gene = "HLA", type = "PROT", suppressMessages = TRUE)
  
  # Find the index of the first sequence that contains the allele name
  match_idx <- grep(allele_name, names(hla_sequences), fixed = TRUE)[1]
  
  # Check if a match was found
  if (is.na(match_idx)) {
    stop("Allele '", allele_name, "' not found in the IMGT/HLA database.")
  }
  
  # Return the sequence as a character string
  as.character(hla_sequences[[match_idx]])
}
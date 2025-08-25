#' Quantify Amino Acid Mismatches Between Two Sequences
#'
#' @description
#' This function compares two amino acid sequences of the same length and
#' returns the number of positions at which the amino acids differ.
#'
#' @param sequence1 A character string representing the first amino acid sequence.
#' @param sequence2 A character string representing the second amino acid sequence.
#'
#' @return An integer representing the number of mismatches.
#'
#' @examples
#' seq1 <- "YFAMYGEKVAHTHVDTLYVRYHY"
#' seq2 <- "YFDMYGEKVAHTHVDTLYVRFHY"
#' quantifyMismatch(seq1, seq2)
#' #> [1] 2
#'
#' @export
quantifyMismatch <- function(sequence1, sequence2) {
  # Input validation
  if (!is.character(sequence1) || !is.character(sequence2)) {
    stop("Input sequences must be character strings.")
  }
  if (nchar(sequence1) != nchar(sequence2)) {
    stop("Input sequences must be of the same length.")
  }

  # Split strings into character vectors
  s1 <- strsplit(sequence1, "")[[1]]
  s2 <- strsplit(sequence2, "")[[1]]

  # Count mismatches
  sum(s1 != s2)
}

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

#' Quantify Eplet Mismatches Between Two Alleles
#'
#' @description
#' This function calculates the number of eplet mismatches between two HLA
#' alleles. It uses the internal `deepMatchR_eplets` dataset.
#'
#' @param allele1 A character string for the first HLA allele (e.g., "A*01:01").
#' @param allele2 A character string for the second HLA allele (e.g., "A*02:01").
#'
#' @return An integer representing the number of mismatched eplets.
#'
#' @importFrom utils data
#' @export
quantifyEpletMismatch <- function(allele1, allele2) {
  # Load the eplet data
  utils::data(deepMatchR_eplets, envir = environment())

  # Convert to data.table and set key for fast subsetting
  eplet_dt <- data.table::as.data.table(deepMatchR_eplets)
  data.table::setkey(eplet_dt, allele)

  # Get eplets for each allele using fast data.table subsetting
  eplets1 <- eplet_dt[.(allele1), eplet, nomatch = 0]
  eplets2 <- eplet_dt[.(allele2), eplet, nomatch = 0]

  # Find the symmetric difference
  mismatched_eplets <- union(setdiff(eplets1, eplets2), setdiff(eplets2, eplets1))

  # Return the count of mismatched eplets
  length(mismatched_eplets)
}

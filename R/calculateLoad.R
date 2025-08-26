#' Calculate Eplet Load Between Donor and Recipient Genotypes
#'
#' @description
#' This function calculates the total eplet load between a donor and a recipient.
#' It iterates over all specified HLA loci, comparing each recipient allele to each
#' donor allele in an exhaustive manner.
#'
#' @param recipient_geno An `hla_genotype` object for the recipient.
#' @param donor_geno An `hla_genotype` object for the donor.
#' @param evidence_level A character vector specifying the evidence levels for
#'   eplets to be included in the analysis. Defaults to `c("A1", "A2")`.
#'
#' @return An integer representing the total eplet load.
#'
#' @examples
#' # Create dummy recipient and donor genotypes
#' recipient <- data.frame(
#'   A_1 = "A*01:01", A_2 = "A*02:01",
#'   B_1 = "B*07:02", B_2 = "B*08:01"
#' )
#' donor <- data.frame(
#'   A_1 = "A*03:01", A_2 = "A*24:02",
#'   B_1 = "B*44:02", B_2 = "B*51:01"
#' )
#'
#' # Convert to hla_genotype objects
#' recipient_geno <- hlaGeno(recipient)
#' donor_geno <- hlaGeno(donor)
#'
#' # Calculate eplet load
#' calculateEpletLoad(recipient_geno, donor_geno)
#'
#' @export
calculateEpletLoad <- function(recipient_geno, donor_geno, evidence_level = c("A1", "A2")) {
  # Input validation
  validateHlaGeno(recipient_geno)
  validateHlaGeno(donor_geno)

  # Get shared loci
  shared_loci <- intersect(recipient_geno$locus_present, donor_geno$locus_present)
  if (length(shared_loci) == 0) {
    stop("No shared loci between donor and recipient.")
  }

  # Get all recipient and donor alleles
  recipient_alleles <- unlist(recipient_geno$data[1, ])
  donor_alleles <- unlist(donor_geno$data[1, ])

  # Remove NA or empty strings
  recipient_alleles <- recipient_alleles[!is.na(recipient_alleles) & recipient_alleles != ""]
  donor_alleles <- donor_alleles[!is.na(donor_alleles) & donor_alleles != ""]

  # Get eplets for all alleles
  utils::data(deepMatchR_eplets, envir = environment())
  eplet_dt <- data.table::as.data.table(deepMatchR_eplets)
  data.table::setkey(eplet_dt, allele)

  if (!is.null(evidence_level)) {
    recipient_eplets <- unique(eplet_dt[allele %in% recipient_alleles & evidence %in% evidence_level, eplet])
    donor_eplets <- unique(eplet_dt[allele %in% donor_alleles & evidence %in% evidence_level, eplet])
  } else {
    recipient_eplets <- unique(eplet_dt[.(recipient_alleles), eplet, nomatch = 0])
    donor_eplets <- unique(eplet_dt[.(donor_alleles), eplet, nomatch = 0])
  }

  # Identify mismatched eplets (present in donor but not in recipient)
  mismatched_eplets <- setdiff(donor_eplets, recipient_eplets)

  return(length(mismatched_eplets))
}

#' Calculate Mismatch Load Between Donor and Recipient Genotypes
#'
#' @description
#' This function calculates the total amino acid mismatch load between a donor
#' and a recipient. It iterates over all specified HLA loci, comparing each
#' recipient allele to each donor allele in an exhaustive manner.
#'
#' @param recipient_geno An `hla_genotype` object for the recipient.
#' @param donor_geno An `hla_genotype` object for the donor.
#'
#' @return An integer representing the total mismatch load.
#'
#' @examples
#' # Create dummy recipient and donor genotypes
#' recipient <- data.frame(
#'   A_1 = "A*01:01", A_2 = "A*02:01",
#'   B_1 = "B*07:02", B_2 = "B*08:01"
#' )
#' donor <- data.frame(
#'   A_1 = "A*01:01", A_2 = "A*03:01",
#'   B_1 = "B*44:02", B_2 = "B*51:01"
#' )
#'
#' # Convert to hla_genotype objects
#' recipient_geno <- hlaGeno(recipient)
#' donor_geno <- hlaGeno(donor)
#'
#' # Calculate mismatch load
#' calculateMismatchLoad(recipient_geno, donor_geno)
#'
#' @export
calculateMismatchLoad <- function(recipient_geno, donor_geno) {
  # Input validation
  validateHlaGeno(recipient_geno)
  validateHlaGeno(donor_geno)

  # Get shared loci
  shared_loci <- intersect(recipient_geno$locus_present, donor_geno$locus_present)
  if (length(shared_loci) == 0) {
    stop("No shared loci between donor and recipient.")
  }

  total_mismatches <- 0

  # Get all recipient and donor alleles
  recipient_alleles <- unlist(recipient_geno$data[1, ])
  donor_alleles <- unlist(donor_geno$data[1, ])

  # Remove NA or empty strings
  recipient_alleles <- recipient_alleles[!is.na(recipient_alleles) & recipient_alleles != ""]
  donor_alleles <- donor_alleles[!is.na(donor_alleles) & donor_alleles != ""]

  # Pre-fetch all unique allele sequences
  all_alleles <- unique(c(recipient_alleles, donor_alleles))
  all_sequences <- sapply(all_alleles, getAlleleSequence)

  # Iterate over each locus
  for (locus in shared_loci) {
    recipient_locus_alleles <- recipient_alleles[grep(paste0("^", locus, "_"), names(recipient_alleles))]
    donor_locus_alleles <- donor_alleles[grep(paste0("^", locus, "_"), names(donor_alleles))]

    # Exhaustive comparison
    for (rec_allele in recipient_locus_alleles) {
      for (don_allele in donor_locus_alleles) {
        if (rec_allele != don_allele) {
          rec_seq <- all_sequences[rec_allele]
          don_seq <- all_sequences[don_allele]
          total_mismatches <- total_mismatches + quantifyMismatch(rec_seq, don_seq)
        }
      }
    }
  }

  return(total_mismatches)
}

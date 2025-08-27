#' Quantify Amino Acid Mismatches With Charge/Polarity Awareness (base R)
#'
#' @description
#' Compares two amino acid sequences (same length) and quantifies mismatches.
#' Each mismatch is annotated for whether it changes residue charge and/or
#' polarity. Users can filter which mismatches to count based on these properties.
#'
#' @param sequence1,sequence2 Character strings of equal length (AAs).
#' @param filter_charge NULL/TRUE/FALSE
#'   - `NULL` (default): ignore charge when filtering (i.e., do not filter).
#'   - `TRUE`: count only mismatches that change charge.
#'   - `FALSE`: count only mismatches that do *not* change charge.
#' @param filter_polarity NULL/TRUE/FALSE.
#'   - `NULL` (default): ignore polarity when filtering.
#'   - `TRUE`: count only mismatches that change polarity.
#'   - `FALSE`: count only mismatches that do *not* change polarity.
#' @param return What to return: one of "count" (default) or "detail"
#' @param na_action One of "exclude" (default), "error", "count". Controls handling
#'   of unknown residues (e.g., X, *, -).
#' @examples
#' seq1 <- "YFAMYGEKVAHTHVDTLYVRYHY"
#' seq2 <- "YFDMYGEKVAHTHVDTLYVRFHY"
#' 
#' # Numerical quantification of mismatches
#' quantifyMismatch(seq1, seq2)
#'
#' # Count only mismatches that change charge
#' quantifyMismatch(seq1, seq2, filter_charge = TRUE)
#'
#' # Count only mismatches that change polarity
#' quantifyMismatch(seq1, seq2, filter_polarity = TRUE)
#'
#' # Count mismatches that change charge AND polarity
#' quantifyMismatch(seq1, seq2, filter_charge = TRUE, filter_polarity = TRUE)
#'
#' @return
#' - If return = "count": integer (count after filters).
#' - If return = "detail": a table with columns:
#'   position, ref, alt, is_mismatch, charge_ref, charge_alt, charge_change,
#'   polarity_ref, polarity_alt, polarity_change, counted.
#'
#' @export
quantifyMismatch <- function(sequence1, 
                             sequence2,
                             filter_charge   = NULL,
                             filter_polarity = NULL,
                             return          = c("count", "detail"),
                             na_action       = c("exclude", "error", "count")) {
  return   <- match.arg(return)
  na_action <- match.arg(na_action)
  
  # --- validation ---
  if (!is.character(sequence1) || !is.character(sequence2)) {
    stop("Input sequences must be character strings.")
  }
  if (length(sequence1) != 1L || length(sequence2) != 1L) {
    stop("Provide single character strings for sequence1 and sequence2.")
  }
  if (nchar(sequence1) != nchar(sequence2)) {
    stop("Input sequences must be of the same length.")
  }
  if (!is.null(filter_charge) && !is.logical(filter_charge)) {
    stop("filter_charge must be NULL, TRUE, or FALSE.")
  }
  if (!is.null(filter_polarity) && !is.logical(filter_polarity)) {
    stop("filter_polarity must be NULL, TRUE, or FALSE.")
  }
  
  # --- AA property maps (physiologic pH) ---
  aa_charge <- c(
    K="pos", R="pos", H="pos",
    D="neg", E="neg",
    A="neu", V="neu", L="neu", I="neu", P="neu",
    M="neu", F="neu", W="neu", G="neu",
    S="neu", T="neu", N="neu", Q="neu", Y="neu", C="neu"
  )
  aa_polarity <- c(
    A="nonpolar", V="nonpolar", L="nonpolar", I="nonpolar", P="nonpolar",
    M="nonpolar", F="nonpolar", W="nonpolar", G="nonpolar",
    S="polar", T="polar", N="polar", Q="polar", Y="polar",
    C="polar", H="polar", K="polar", R="polar", D="polar", E="polar"
  )
  
  # --- split & normalize ---
  s1 <- strsplit(toupper(sequence1), "", fixed = TRUE)[[1]]
  s2 <- strsplit(toupper(sequence2), "", fixed = TRUE)[[1]]
  
  standard <- names(aa_charge)
  s1_known <- s1 %in% standard
  s2_known <- s2 %in% standard
  any_unknown <- any(!s1_known | !s2_known)
  
  if (any_unknown && na_action == "error") {
    bad_pos <- which(!s1_known | !s2_known)
    stop(sprintf("Unknown/unsupported residue(s) at position(s): %s",
                 paste(bad_pos, collapse = ", ")))
  }
  
  # --- mismatch flags & properties ---
  is_mismatch <- s1 != s2
  
  charge_ref <- unname(aa_charge[s1]); charge_ref[is.na(charge_ref)] <- NA
  charge_alt <- unname(aa_charge[s2]); charge_alt[is.na(charge_alt)] <- NA
  pol_ref    <- unname(aa_polarity[s1]); pol_ref[is.na(pol_ref)] <- NA
  pol_alt    <- unname(aa_polarity[s2]); pol_alt[is.na(pol_alt)] <- NA
  
  charge_change   <- ifelse(is_mismatch, charge_ref != charge_alt, FALSE)
  polarity_change <- ifelse(is_mismatch, pol_ref    != pol_alt,    FALSE)
  
  # Handle unknowns according to na_action
  if (any_unknown) {
    unk <- (!s1_known | !s2_known)
    if (na_action == "exclude") {
      charge_change[unk]   <- NA
      polarity_change[unk] <- NA
    } else if (na_action == "count") {
      # keep mismatch counting; property deltas remain NA where unknown
      charge_change[unk & is_mismatch]   <- NA
      polarity_change[unk & is_mismatch] <- NA
    }
  }
  
  # --- apply filters ---
  counted <- is_mismatch
  if (!is.null(filter_charge)) {
    counted <- counted & (charge_change %in% filter_charge)
  }
  if (!is.null(filter_polarity)) {
    counted <- counted & (polarity_change %in% filter_polarity)
  }
  if (na_action == "exclude") {
    if (!is.null(filter_charge))   counted <- counted & !is.na(charge_change)
    if (!is.null(filter_polarity)) counted <- counted & !is.na(polarity_change)
  }
  
  # Fast path: return just the count
  if (return == "count") {
    return(sum(counted, na.rm = TRUE))
  }
  
  # Build a base data.frame (no extra deps)
  detail <- data.frame(
    position        = seq_along(s1),
    ref             = s1,
    alt             = s2,
    is_mismatch     = is_mismatch,
    charge_ref      = charge_ref,
    charge_alt      = charge_alt,
    charge_change   = as.logical(charge_change),
    polarity_ref    = pol_ref,
    polarity_alt    = pol_alt,
    polarity_change = as.logical(polarity_change),
    counted         = as.logical(counted),
    stringsAsFactors = FALSE
  )
  
  if (return == "detail") {
    return(detail)
  }
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
#' @param evidence_level Character vector indicating the antibody reactivity levels to keep.
#'   Defaults to \code{c("A1", "A2")}, which represent antibody-confirmed eplets.
#'   Other acceptable levels include \code{"B"}, \code{"D"}, or \code{NULL} to apply no filter.
#'
#' @return An integer representing the number of mismatched eplets.
#'
#' @importFrom utils data
#' @importFrom data.table as.data.table setkey
#' @export
quantifyEpletMismatch <- function(allele1,
                                  allele2,
                                  evidence_level = c("A1", "A2")) {
  # Load the eplet data
  utils::data(deepMatchR_eplets, envir = environment())

  # Convert to data.table and set key for fast subsetting
  eplet_dt <- data.table::as.data.table(deepMatchR_eplets)
  data.table::setkey(eplet_dt, allele)


  
  if(!is.null(evidence_level)) {
    eplets1 <- eplet_dt[
           allele %in% allele1 & evidence %in% evidence_level,
           eplet]
    eplets2 <- eplet_dt[
      allele %in% allele2 & evidence %in% evidence_level,
      eplet]
  } else {
    # Get eplets for each allele using fast data.table subsetting
    eplets1 <- eplet_dt[.(allele1), eplet, nomatch = 0]
    eplets2 <- eplet_dt[.(allele2), eplet, nomatch = 0]
  }
  
 

  # Find the symmetric difference
  mismatched_eplets <- union(setdiff(eplets1, eplets2), setdiff(eplets2, eplets1))

  # Return the count of mismatched eplets
  length(mismatched_eplets)
}

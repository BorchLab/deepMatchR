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
#' @param type Character string. The type of alignment to perform. Defaults to 
#'   `"global"` but allows for `"local"` and `"overlap"`-based alignments of the 
#'   sequences.
#' @param substitutionMatrix Character string or numeric matrix. Substitution 
#'   scoring matrix used during sequence alignment. Defaults to `"BLOSUM80"`, 
#'   which provides a conservative amino acid similarity scheme suitable for 
#'   closely related protein sequences.
#' @param gapOpening Numeric scalar. Penalty score applied when initiating a 
#'   new gap during alignment. Higher values discourage insertion/deletion 
#'   events and yield more contiguous alignments. Default is `10`.
#' @param gapExtension Numeric scalar. Penalty score applied when extending an 
#'   existing gap. Smaller values permit longer continuous gaps, while larger 
#'   values favor shorter gaps. Default is `1`.
#' @param count_gaps Logical (default `TRUE`). If `TRUE`, positions where one 
#'   sequence contains a gap (`-`) and the other contains an amino acid are 
#'   treated as mismatches and included in the mismatch count. 
#' 
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
#' @importFrom pwalign pattern subject pairwiseAlignment score
#' @export
quantifyMismatch <- function(sequence1, 
                                     sequence2,
                                     filter_charge   = NULL,
                                     filter_polarity = NULL,
                                     return          = c("count", "detail"),
                                     na_action       = c("exclude", "error", "count"),
                                     type            = c("global", "local", "overlap"),
                                     substitutionMatrix = "BLOSUM80",
                                     gapOpening      = 10,
                                     gapExtension    = 1,
                                     count_gaps      = TRUE) {
  return    <- match.arg(return)
  na_action <- match.arg(na_action)
  type      <- match.arg(type)
  
  if (!is.character(sequence1) || !is.character(sequence2) ||
      length(sequence1) != 1L || length(sequence2) != 1L) {
    stop("Provide single character strings for sequence1 and sequence2.")
  }
  if (!is.null(filter_charge) && !is.logical(filter_charge)) {
    stop("filter_charge must be NULL, TRUE, or FALSE.")
  }
  if (!is.null(filter_polarity) && !is.logical(filter_polarity)) {
    stop("filter_polarity must be NULL, TRUE, or FALSE.")
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings", " not installed, install or choose a different `method`.",
         call. = FALSE)
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
  standard <- names(aa_charge)
  
  # --- Align (global by default) ---
  data(list = substitutionMatrix, package = "pwalign", envir = environment())
  mat <- if (is.character(substitutionMatrix)) get(substitutionMatrix) else substitutionMatrix
  
  pwa <- pwalign::pairwiseAlignment(
    pattern = Biostrings::AAString(sequence1),
    subject = Biostrings::AAString(sequence2),
    type = type,
    substitutionMatrix = mat,
    gapOpening   = gapOpening,
    gapExtension = gapExtension
  )
  
  aln1 <- as.character(pwalign::pattern(pwa))   # with '-'
  aln2 <- as.character(pwalign::subject(pwa))   # with '-'
  
  s1 <- strsplit(aln1, "", fixed = TRUE)[[1]]
  s2 <- strsplit(aln2, "", fixed = TRUE)[[1]]
  
  # --- unknowns & gaps ---
  is_gap1 <- s1 == "-"
  is_gap2 <- s2 == "-"
  is_gap  <- is_gap1 | is_gap2
  
  s1_known <- s1 %in% standard
  s2_known <- s2 %in% standard
  any_unknown <- any((!s1_known & !is_gap1) | (!s2_known & !is_gap2))
  
  if (any_unknown && na_action == "error") {
    bad_pos <- which((!s1_known & !is_gap1) | (!s2_known & !is_gap2))
    stop(sprintf("Unknown/unsupported residue(s) at aligned position(s): %s",
                 paste(bad_pos, collapse = ", ")))
  }
  
  # --- mismatches (subs vs gaps) ---
  is_mismatch <- s1 != s2
  # If we choose to ignore gaps for counting:
  if (!count_gaps) {
    is_mismatch[is_gap] <- FALSE
  }
  
  # --- charge / polarity annotations (NA for gaps) ---
  charge_ref <- ifelse(is_gap1, NA, unname(aa_charge[s1]))
  charge_alt <- ifelse(is_gap2, NA, unname(aa_charge[s2]))
  pol_ref    <- ifelse(is_gap1, NA, unname(aa_polarity[s1]))
  pol_alt    <- ifelse(is_gap2, NA, unname(aa_polarity[s2]))
  
  charge_change   <- ifelse(is_mismatch & !is_gap, charge_ref != charge_alt, FALSE)
  polarity_change <- ifelse(is_mismatch & !is_gap, pol_ref    != pol_alt,    FALSE)
  
  # Handle unknowns per na_action
  if (any_unknown) {
    unk <- ((!s1_known & !is_gap1) | (!s2_known & !is_gap2))
    if (na_action == "exclude") {
      charge_change[unk]   <- NA
      polarity_change[unk] <- NA
    } else if (na_action == "count") {
      # keep mismatch counting; deltas remain NA where unknown
      charge_change[unk & is_mismatch]   <- NA
      polarity_change[unk & is_mismatch] <- NA
    }
  }
  
  # --- apply filters ---
  counted <- is_mismatch
  if (!is.null(filter_charge)) {
    # charge filter applies only to non-gap substitutions; gaps have NA charge_change
    counted <- counted & (is_gap | (charge_change %in% filter_charge))
  }
  if (!is.null(filter_polarity)) {
    counted <- counted & (is_gap | (polarity_change %in% filter_polarity))
  }
  if (na_action == "exclude") {
    if (!is.null(filter_charge))   counted <- counted & (is_gap | !is.na(charge_change))
    if (!is.null(filter_polarity)) counted <- counted & (is_gap | !is.na(polarity_change))
  }
  
  if (return == "count") {
    return(sum(counted, na.rm = TRUE))
  }
  
  # --- detail table ---
  detail <- data.frame(
    alignment_position = seq_along(s1),
    ref                = s1,
    alt                = s2,
    is_gap_ref         = is_gap1,
    is_gap_alt         = is_gap2,
    is_mismatch        = is_mismatch,
    charge_ref         = charge_ref,
    charge_alt         = charge_alt,
    charge_change      = as.logical(charge_change),
    polarity_ref       = pol_ref,
    polarity_alt       = pol_alt,
    polarity_change    = as.logical(polarity_change),
    counted            = as.logical(counted),
    stringsAsFactors   = FALSE
  )
  
  attr(detail, "alignment") <- list(
    aligned_ref = aln1,
    aligned_alt = aln2,
    score       = pwalign::score(pwa),
    type        = type
  )
  
  detail
}

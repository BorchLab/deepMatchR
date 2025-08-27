#' Calculate Mismatch Load Between Donor and Recipient Genotypes
#'
#' @description
#' Computes amino-acid mismatch burden between donor and recipient across HLA loci
#' by calling `quantifyMismatch()` on each recipient–donor allele pair per locus.
#' Can return (1) a single total, (2) a per-locus summary table, or
#' (3) a pairwise recipient-vs-donor allele matrix for one locus.
#'
#' @param recipient_geno An `hla_genotype` object for the recipient.
#' @param donor_geno An `hla_genotype` object for the donor.
#' @param loci Character vector of loci to include (e.g., c("A","B","C")), or NULL
#'   to use all shared loci.
#' @param filter_charge NULL/TRUE/FALSE. Passed to [quantifyMismatch()].
#' @param filter_polarity NULL/TRUE/FALSE. Passed to [quantifyMismatch()].
#' @param na_action One of "exclude" (default), "error", "count". Passed through.
#' @param return What to return: "total" (default), "per_locus", or "pairwise".
#' @param pairwise_locus When `return = "pairwise"`, the locus to visualize
#'   (e.g., "A", "B", "C"). Ignored otherwise.
#'
#' @return
#' - If `return = "total"`: single integer total mismatch load.
#' - If `return = "per_locus"`: data.frame with columns `locus`, `mismatch_load`.
#' - If `return = "pairwise"`: numeric matrix of pairwise counts (rows = recipient alleles,
#'   columns = donor alleles) for `pairwise_locus`.
#'
#' @examples
#' # Toy genotypes
#' recipient <- data.frame(
#'   A_1 = "A*01:01", A_2 = "A*02:01",
#'   B_1 = "B*07:02", B_2 = "B*08:01"
#' )
#' donor <- data.frame(
#'   A_1 = "A*01:01", A_2 = "A*03:01",
#'   B_1 = "B*44:02", B_2 = "B*51:01"
#' )
#' rgeno <- hlaGeno(recipient)
#' dgeno <- hlaGeno(donor)
#'
#' # 1) Total load (all loci, all mismatches)
#' calculateMismatchLoad(rgeno, dgeno)
#'
#' # 2) Per-locus summary
#' calculateMismatchLoad(rgeno, dgeno, return = "per_locus")
#'
#' # 3) Pairwise matrix for B locus
#' mB <- calculateMismatchLoad(rgeno, dgeno, return = "pairwise", pairwise_locus = "B")
#' mB
#'
#' # 4) Apply biophysical filters
#' calculateMismatchLoad(rgeno, dgeno, filter_charge = TRUE)                # charge-changing only
#' calculateMismatchLoad(rgeno, dgeno, filter_polarity = TRUE)              # polarity-changing only
#' calculateMismatchLoad(rgeno, dgeno, filter_charge = TRUE, filter_polarity = TRUE)  # both
#'
#' @export
calculateMismatchLoad <- function(recipient_geno,
                                  donor_geno,
                                  loci = NULL,
                                  filter_charge = NULL,
                                  filter_polarity = NULL,
                                  na_action = c("exclude", "error", "count"),
                                  return = c("total", "per_locus", "pairwise"),
                                  pairwise_locus = NULL) {
  # --- validation ---
  na_action <- match.arg(na_action)
  return    <- match.arg(return)
  
  validateHlaGeno(recipient_geno)
  validateHlaGeno(donor_geno)
  
  shared_loci <- intersect(recipient_geno$locus_present, donor_geno$locus_present)
  if (!is.null(loci)) shared_loci <- intersect(shared_loci, loci)
  if (length(shared_loci) == 0) stop("No shared loci between donor and recipient (after filtering 'loci').")
  
  if (return == "pairwise") {
    if (length(pairwise_locus) != 1L || !is.character(pairwise_locus)) {
      stop("When return = 'pairwise', provide pairwise_locus as a single character locus (e.g., 'B').")
    }
    if (!(pairwise_locus %in% shared_loci)) {
      stop(sprintf("pairwise_locus '%s' is not in shared_loci: %s",
                   pairwise_locus, paste(shared_loci, collapse = ", ")))
    }
  }
  
  # --- collect allele labels (1st row of $data), optionally subset by loci pattern ---
  get_idx_for_loci <- function(g, loci_vec) {
    if (is.null(loci_vec)) return(seq_len(ncol(g$data)))
    # Columns typically like "A_1", "B_2"; select those starting with e.g. "A_"
    pat <- paste0("^(", paste(loci_vec, collapse = "|"), ")_")
    grep(pat, colnames(g$data))
  }
  
  r_idx <- get_idx_for_loci(recipient_geno, shared_loci)
  d_idx <- get_idx_for_loci(donor_geno,     shared_loci)
  
  r_alleles <- unlist(recipient_geno$data[1, r_idx, drop = FALSE])
  d_alleles <- unlist(donor_geno$data[1,     d_idx, drop = FALSE])
  
  # Remove NA/empty
  r_alleles <- r_alleles[!is.na(r_alleles) & nzchar(r_alleles)]
  d_alleles <- d_alleles[!is.na(d_alleles) & nzchar(d_alleles)]
  
  # --- Prefetch sequences once for speed ---
  all_alleles <- unique(c(r_alleles, d_alleles))
  if (length(all_alleles) == 0) stop("No allele strings found in genotype data.")
  seq_map <- sapply(all_alleles, getAlleleSequence, USE.NAMES = TRUE)
  
  # Helper: compute total for one locus
  locus_total <- function(L) {
    rL <- r_alleles[grep(paste0("^", L, "_"), names(r_alleles))]
    dL <- d_alleles[grep(paste0("^", L, "_"), names(d_alleles))]
    if (length(rL) == 0L || length(dL) == 0L) return(0L)
    
    total <- 0L
    for (ra in rL) for (da in dL) {
      if (identical(ra, da)) next
      total <- total + quantifyMismatch(
        sequence1 = unname(seq_map[ra]),
        sequence2 = unname(seq_map[da]),
        filter_charge   = filter_charge,
        filter_polarity = filter_polarity,
        na_action       = na_action,
        return          = "count"
      )
    }
    total
  }
  
  if (return == "total") {
    # Sum over loci
    loads <- vapply(shared_loci, locus_total, FUN.VALUE = integer(1))
    return(as.integer(sum(loads)))
  }
  
  if (return == "per_locus") {
    loads <- vapply(shared_loci, locus_total, FUN.VALUE = integer(1))
    return(data.frame(locus = shared_loci,
                      mismatch_load = as.integer(loads),
                      row.names = NULL, check.names = FALSE))
  }
  
  # return == "pairwise"
  L <- pairwise_locus
  rL <- r_alleles[grep(paste0("^", L, "_"), names(r_alleles))]
  dL <- d_alleles[grep(paste0("^", L, "_"), names(d_alleles))]
  if (length(rL) == 0L || length(dL) == 0L) stop("No alleles found for locus ", L)
  
  mat <- outer(
    rL, dL,
    Vectorize(function(ra, da) {
      if (identical(ra, da)) return(0L)
      quantifyMismatch(
        sequence1 = unname(seq_map[ra]),
        sequence2 = unname(seq_map[da]),
        filter_charge   = filter_charge,
        filter_polarity = filter_polarity,
        na_action       = na_action,
        return          = "count"
      )
    })
  )
  dimnames(mat) <- list(recipient = rL, donor = dL)
  return(mat)
}

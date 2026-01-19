#' Calculate Eplet Load Between Donor and Recipient Genotypes
#'
#' @description
#' Aggregates **donor-specific eplets** (present in donor, absent in recipient)
#' over loci. Can return a single total, a per-locus summary, or a pairwise
#' recipient-vs-donor allele matrix for a chosen locus. Supports filtering by
#' **evidence**, **exposition**, and **reactivity**.
#'
#' @param recipient_geno,donor_geno `hla_genotype` objects.
#' @param loci Character vector of loci to include (e.g. `c("A","B","C")`), or
#'   `NULL` (default) for all shared loci.
#' @param evidence_level Character vector of evidence levels to include
#'   (default `c("A1","A2")`). Use `NULL` for no evidence filter.
#' @param exposition_filter Character vector of exposition classes to include;
#'   `NULL` for no filter.
#' @param reactivity_filter Character vector of reactivity classes to include;
#'   `NULL` for no filter.
#' @param return What to return: `"total"` (default), `"per_locus"`, or `"pairwise"`.
#' @param pairwise_locus When `return = "pairwise"`, the single locus to compute
#'   (e.g., `"B"`). Ignored otherwise.
#'
#' @return
#' - If `return="total"`: integer total donor-specific eplet load.
#' - If `return="per_locus"`: `data.frame` with columns `locus`, `eplet_load`.
#' - If `return="pairwise"`: numeric matrix where entry `[i,j]` is the count of
#'   donor-specific eplets for donor allele `j` vs recipient allele `i` at `pairwise_locus`.
#'
#' @examples
#' # Dummy genotypes
#' recipient <- data.frame(
#'   A_1 = "A*01:01", A_2 = "A*02:01",
#'   B_1 = "B*07:02", B_2 = "B*08:01"
#' )
#' donor <- data.frame(
#'   A_1 = "A*03:01", A_2 = "A*24:02",
#'   B_1 = "B*44:02", B_2 = "B*51:01"
#' )
#' rgeno <- hlaGeno(recipient); dgeno <- hlaGeno(donor)
#'
#' # Total donor-specific eplet load (A1/A2 evidence)
#' calculateEpletLoad(rgeno, dgeno)
#'
#' # Per locus
#' calculateEpletLoad(rgeno, dgeno, return = "per_locus")
#'
#' # Pairwise matrix for B locus (rows=recipient alleles, cols=donor alleles)
#' mB <- calculateEpletLoad(rgeno, dgeno, return = "pairwise", pairwise_locus = "B")
#' mB
#'
#' # Apply additional filters
#' calculateEpletLoad(rgeno, dgeno, exposition_filter = "High")
#' calculateEpletLoad(rgeno, dgeno, reactivity_filter = c("IgG"))
#'
#' @importFrom data.table as.data.table setkey
#' @export
calculateEpletLoad <- function(recipient_geno,
                               donor_geno,
                               loci = NULL,
                               evidence_level    = c("A1", "A2"),
                               exposition_filter = NULL,
                               reactivity_filter = NULL,
                               return = c("total", "per_locus", "pairwise"),
                               pairwise_locus = NULL) {
  return <- match.arg(return)
  
  # Validate genotypes
  validateHlaGeno(recipient_geno)
  validateHlaGeno(donor_geno)
  
  # Determine loci
  shared_loci <- intersect(recipient_geno$locus_present, donor_geno$locus_present)
  if (!is.null(loci)) shared_loci <- intersect(shared_loci, loci)
  if (length(shared_loci) == 0L) stop("No shared loci between donor and recipient after filtering 'loci'.")
  
  # Extract allele strings (1st row), subset columns by loci
  pick_cols <- function(g, loci_vec) {
    if (is.null(loci_vec)) return(seq_len(ncol(g$data)))
    grep(paste0("^(", paste(loci_vec, collapse = "|"), ")_"), colnames(g$data))
  }
  r_idx <- pick_cols(recipient_geno, shared_loci)
  d_idx <- pick_cols(donor_geno,     shared_loci)
  
  r_alleles <- unlist(recipient_geno$data[1, r_idx, drop = FALSE])
  d_alleles <- unlist(donor_geno$data[1, d_idx, drop = FALSE])
  
  # Clean
  r_alleles <- r_alleles[!is.na(r_alleles) & nzchar(r_alleles)]
  d_alleles <- d_alleles[!is.na(d_alleles) & nzchar(d_alleles)]
  if (length(r_alleles) == 0L || length(d_alleles) == 0L) stop("No allele strings found in genotype data.")
  
  # Load eplet table & apply global filters once
  eplet_dt <- data.table::as.data.table(deepMatchR_eplets)
  
  # Apply conjunctive filters if provided
  if (!is.null(evidence_level))    eplet_dt <- eplet_dt[evidence   %in% evidence_level]
  if (!is.null(exposition_filter)) eplet_dt <- eplet_dt[exposition %in% exposition_filter]
  if (!is.null(reactivity_filter)) eplet_dt <- eplet_dt[reactivity %in% reactivity_filter]
  
  # Index for fast lookup
  data.table::setkey(eplet_dt, allele)
  
  # Precompute eplet sets per allele (after filters)
  all_alleles <- unique(c(r_alleles, d_alleles))
  eplet_sets <- lapply(all_alleles, function(a)
    unique(eplet_dt[list(a), on = "allele", eplet, nomatch = 0L])
  )
  names(eplet_sets) <- all_alleles
  
  # Helper: donor-specific eplet count for one pair
  pair_count <- function(rec_a, don_a) {
    if (identical(rec_a, don_a)) return(0L)
    rec_set <- eplet_sets[[rec_a]]
    don_set <- eplet_sets[[don_a]]
    length(setdiff(don_set, rec_set))
  }
  
  # Helper: total per locus
  locus_total <- function(L) {
    rL <- r_alleles[grep(paste0("^", L, "_"), names(r_alleles))]
    dL <- d_alleles[grep(paste0("^", L, "_"), names(d_alleles))]
    if (length(rL) == 0L || length(dL) == 0L) return(0L)
    total <- 0L
    for (ra in rL) for (da in dL) total <- total + pair_count(ra, da)
    total
  }
  
  if (return == "total") {
    loads <- vapply(shared_loci, locus_total, FUN.VALUE = integer(1))
    return(as.integer(sum(loads)))
  }
  
  if (return == "per_locus") {
    loads <- vapply(shared_loci, locus_total, FUN.VALUE = integer(1))
    return(data.frame(locus = shared_loci,
                      eplet_load = as.integer(loads),
                      row.names = NULL, check.names = FALSE))
  }
  
  # pairwise
  if (length(pairwise_locus) != 1L || !is.character(pairwise_locus))
    stop("When return = 'pairwise', provide pairwise_locus as a single character (e.g., 'B').")
  if (!(pairwise_locus %in% shared_loci))
    stop(sprintf("pairwise_locus '%s' is not shared; shared loci: %s",
                 pairwise_locus, paste(shared_loci, collapse = ", ")))
  
  rL <- r_alleles[grep(paste0("^", pairwise_locus, "_"), names(r_alleles))]
  dL <- d_alleles[grep(paste0("^", pairwise_locus, "_"), names(d_alleles))]
  if (length(rL) == 0L || length(dL) == 0L) stop("No alleles found for locus ", pairwise_locus)
  
  mat <- outer(
    rL, dL,
    Vectorize(function(ra, da) pair_count(ra, da))
  )
  dimnames(mat) <- list(recipient = rL, donor = dL)
  mat
}

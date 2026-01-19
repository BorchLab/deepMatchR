#' Quantify Eplet Mismatches Between Two Alleles
#'
#' @description
#' Counts the number of **eplet mismatches** between two HLA alleles, using the
#' internal `deepMatchR_eplets` dataset. Supports filtering by **evidence**
#' (e.g., A1/A2), **exposition** (e.g., "High", "Moderate"), and **reactivity**
#' (e.g., IgG/IgM labels if present in the table).
#'
#' @param allele1,allele2 Character scalars, e.g. `"A*01:01"`, `"A*02:01"`.
#' @param evidence_level Character vector of evidence levels to include
#'   (default `c("A1","A2")`). Use `NULL` for no evidence filter.
#' @param exposition_filter Character vector of exposition classes to include
#'   (e.g., `c("High","Moderate")`). Default `NULL` = no exposition filter.
#' @param reactivity_filter Character vector of reactivity classes to include
#'   (dataset-dependent). Default `NULL` = no reactivity filter.
#'
#' @return Integer: size of the **symmetric difference** of eplet sets between
#'   the two alleles after filters (i.e., eplets present in one allele but not the other).
#'
#' @examples
#' # Count eplet mismatches between two alleles
#' quantifyEpletMismatch("A*01:01", "A*02:01")
#'
#' # With evidence level filter
#' quantifyEpletMismatch("A*01:01", "A*02:01", evidence_level = "A1")
#'
#' # Same allele returns 0
#' quantifyEpletMismatch("A*01:01", "A*01:01")
#'
#' @importFrom data.table as.data.table setkey
#' @export
quantifyEpletMismatch <- function(allele1,
                                  allele2,
                                  evidence_level   = c("A1", "A2"),
                                  exposition_filter = NULL,
                                  reactivity_filter = NULL) {
  # Fast exits
  if (!is.character(allele1) || !is.character(allele2) ||
      length(allele1) != 1L || length(allele2) != 1L) {
    stop("allele1 and allele2 must be length-1 character strings.")
  }
  if (identical(allele1, allele2)) return(0L)
  
  # Load & prep
  eplet_dt <- data.table::as.data.table(deepMatchR_eplets)
  data.table::setkey(eplet_dt, allele)
  
  # Build conjunctive filter
  filt <- rep_len(TRUE, nrow(eplet_dt))
  if (!is.null(evidence_level))    filt <- filt & (eplet_dt$evidence   %in% evidence_level)
  if (!is.null(exposition_filter)) filt <- filt & (eplet_dt$exposition %in% exposition_filter)
  if (!is.null(reactivity_filter)) filt <- filt & (eplet_dt$reactivity %in% reactivity_filter)
  
  # Subset once
  eplt <- eplet_dt[filt]
  
  # Eplet sets per allele
  e1 <- unique(eplt[list(allele1), on = "allele", eplet, nomatch = 0L])
  e2 <- unique(eplt[list(allele2), on = "allele", eplet, nomatch = 0L])
  
  # Symmetric difference size
  length(union(setdiff(e1, e2), setdiff(e2, e1)))
}
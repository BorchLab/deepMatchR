#' CREG–Allele Mapping Data
#'
#' A curated mapping of IMGT/HLA allele strings to their corresponding
#' serologic antigen assignments and Cross-Reactive Groups (CREGs). This table
#' is intended for convenience functions that summarize antibody specificity,
#' collapse alleles to serology, or group responses by CREG for reporting.
#'
#' @details
#' - **Allele strings** follow IMGT/HLA nomenclature (e.g., `"A*02:01"`).
#' - **Serology** uses conventional two-digit antigen labels (e.g., `"A2"`,
#'   `"B8"`). Depending on your workflow, you may want to harmonize with
#'   vendor-specific naming.
#' - **CREG** indicates the cross-reactive group label used in many transplant
#'   workflows to approximate serologic cross-reactivity.
#'
#' @format A data frame with the following variables:
#' \describe{
#'   \item{allele}{`character`. IMGT/HLA allele (e.g., `"B*07:02"`).}
#'   \item{serology}{`character`. Serologic antigen assignment (e.g., `"B7"`).}
#'   \item{CREG}{`character`. Cross-Reactive Group label (e.g., `"CREG07"`).}
#' }
#'
#' @section Typical use:
#' - Collapsing allele-level reactivity to serology/CREG.
#' - Building summary tables/plots by serology or CREG.
#'
#' @usage data(deepMatchR_cregs)
#' @seealso \code{\link{deepMatchR_eplets}}, \code{\link{deepMatchR_example}}
#' @keywords datasets HLA serology CREG
#' @note Source is an internal curation aligned to common CREG practice; see the
#'   package vignette for curation notes and limitations. Always verify against
#'   your lab’s approved references.
"deepMatchR_cregs"



#' HLA Eplet Assignments (Registry-derived)
#'
#' Per-allele eplet annotations derived from the HLA Epitope Registry, filtered
#' to retain evidence classes A1, A2, B, and D. Each row links an eplet to an
#' HLA allele along with Registry metadata fields commonly used for analysis.
#'
#' @details
#' **Field meanings (as used in this package):**
#' \itemize{
#'   \item \strong{eplet}: Short alphanumeric eplet identifier (e.g., `"82LR"`, `"1C"`).
#'   \item \strong{exposition}: Qualitative description of surface exposure or
#'         structural context (e.g., `"High"`). This reflects the Registry’s
#'         exposition field.
#'   \item \strong{reactivity}: Free-text/flag describing observed reactivity
#'         patterns in the Registry (often `NA` if not specified).
#'   \item \strong{evidence}: Evidence class label from the Registry; records in
#'         this dataset are filtered to A1, A2, B, or D.
#'   \item \strong{allele}: IMGT/HLA allele string to which the eplet is
#'         assigned (e.g., `"B*07:02"`).
#' }
#'
#' @section Caveats:
#' - The Registry is a living resource; re-download or update your cache
#'   regularly for production use.
#' - Evidence codes indicate strength/quality of support and are not equivalent
#'   to clinical validity. Apply your lab’s validation and cutoffs.
#'
#' @format A data frame with 5 variables:
#' \describe{
#'   \item{eplet}{`character`. Eplet identifier.}
#'   \item{exposition}{`character`. Structural/surface exposure categorization.}
#'   \item{reactivity}{`character` or `NA`. Registry reactivity note/flag.}
#'   \item{evidence}{`character`. Evidence class; subset of \{A1, A2, B, D\}.}
#'   \item{allele}{`character`. IMGT/HLA allele string.}
#' }
#'
#' @usage data(deepMatchR_eplets)
#' @source HLA Epitope Registry (\url{https://www.epregistry.com.br/}); processed
#'   and filtered by the package authors for reproducible analyses.
#' @seealso \code{\link{deepMatchR_cregs}}
#' @keywords datasets HLA eplets epitope
"deepMatchR_eplets"



#' Example SAB (Class I/II) and PRA Panels
#'
#' A small, fixed-format example object demonstrating the input structure used by
#' this package’s utilities for single antigen bead (SAB) Class I / Class II and
#' panel reactive antibody (PRA) data. Useful for examples, vignettes, and unit
#' tests without requiring PHI or proprietary vendor exports.
#'
#' @details
#' \strong{Structure:} a named list of length 3:
#' \itemize{
#'   \item \code{ClassI_example}: data frame with columns
#'     \code{BeadID} (`integer`), \code{SpecAbbr} (`character`, antigen abbreviations),
#'     \code{Specificity} (`character`, comma-delimited allele list),
#'     \code{NormalValue} (`numeric`, normalized MFI or vendor-provided normalization),
#'     \code{RawData} (`numeric`, raw MFI), \code{CountValue} (`integer`, bead/event count).
#'   \item \code{ClassII_example}: same schema as Class I, but for class II
#'     specificities (e.g., DR, DQ, DP).
#'   \item \code{PRA}: data frame with columns
#'     \code{BeadID} (`integer`), \code{SpecAbbr} (`character`),
#'     \code{Specificity} (`character`), \code{NormalValue} (`numeric` or `NA`),
#'     \code{RawData} (`numeric`), \code{CountValue} (`integer`).
#' }
#'
#' \strong{Notes:}
#' - Column names are kept vendor-agnostic but mimic common exports.
#' - \code{SpecAbbr} values are comma-separated antigen abbreviations with
#'   padding dashes, e.g., `"A2,-,-,-,..."`.
#' - \code{Specificity} values are comma-separated IMGT/HLA alleles aligned to
#'   the abbreviations in \code{SpecAbbr}.
#' - \code{NormalValue} may be \code{NA} when not supplied by the instrument
#'   export or when illustrative only.
#'
#' @format A named list of length 3 containing data frames as described above.
#'
#' @usage data(deepMatchR_example)
#' @seealso \code{\link{deepMatchR_eplets}}, \code{\link{deepMatchR_cregs}}
#' @keywords datasets SAB PRA MFI
"deepMatchR_example"

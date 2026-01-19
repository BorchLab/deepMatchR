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


#' WMDA DNA-to-Serology Mapping
#'
#' A lookup table mapping HLA alleles to their serological equivalents based on
#' WMDA (World Marrow Donor Association) nomenclature. Used internally by
#' \code{\link{toSerology}} for allele-to-serology conversion.
#'
#' @format A data.table with the following columns:
#' \describe{
#'   \item{locus}{`character`. HLA locus with asterisk (e.g., `"A*"`, `"B*"`, `"DRB1*"`).}
#'   \item{allele_2f}{`character`. Two-field allele designation (e.g., `"01:01"`, `"07:02"`).}
#'   \item{serology}{`character`. Serological antigen number (e.g., `"1"`, `"7"`).}
#' }
#'
#' @details
#' The mapping prioritizes serology assignments in this order:
#' \enumerate{
#'   \item Unambiguous assignments
#'   \item Possible assignments
#'   \item Assumed assignments
#'   \item Expert assignments
#' }
#'
#' @source WMDA nomenclature files from IMGT/HLA GitHub repository
#'   (\url{https://github.com/ANHIG/IMGTHLA/tree/Latest/wmda}).
#'
#' @usage data(deepMatchR_wmda_serology)
#' @seealso \code{\link{toSerology}}, \code{\link{updateWmdaData}},
#'   \code{\link{deepMatchR_wmda_splits}}, \code{\link{deepMatchR_wmda_pgroups}}
#' @keywords datasets HLA serology WMDA
"deepMatchR_wmda_serology"


#' WMDA Broad-to-Split Antigen Relationships
#'
#' A lookup table mapping broad serological antigens to their split antigens
#' based on WMDA nomenclature. Used internally by \code{\link{toSerology}}
#' when \code{resolve_splits = TRUE}.
#'
#' @format A data.table with the following columns:
#' \describe{
#'   \item{locus}{`character`. Serology locus prefix (e.g., `"A"`, `"B"`, `"DR"`).}
#'   \item{broad}{`character`. Broad antigen number (e.g., `"2"`, `"5"`).}
#'   \item{splits}{`character`. Pipe-separated split antigen numbers (e.g., `"15|16"`).}
#' }
#'
#' @details
#' Common examples of broad-to-split relationships:
#' \itemize{
#'   \item DR2 -> DR15, DR16
#'   \item DR5 -> DR11, DR12
#'   \item B5 -> B51, B52
#' }
#'
#' @source WMDA nomenclature files from IMGT/HLA GitHub repository
#'   (\url{https://github.com/ANHIG/IMGTHLA/tree/Latest/wmda}).
#'
#' @usage data(deepMatchR_wmda_splits)
#' @seealso \code{\link{toSerology}}, \code{\link{deepMatchR_wmda_serology}}
#' @keywords datasets HLA serology WMDA
"deepMatchR_wmda_splits"


#' WMDA P-Group Definitions
#'
#' A lookup table containing HLA P-group definitions from WMDA nomenclature.
#' P-groups are sets of alleles with identical protein sequences in the
#' antigen recognition site. Used internally by \code{\link{toSerology}}
#' to resolve P-group notation.
#'
#' @format A data.table with the following columns:
#' \describe{
#'   \item{locus}{`character`. HLA locus (e.g., `"A"`, `"B"`, `"DRB1"`).}
#'   \item{p_group}{`character`. P-group designation (e.g., `"01:01P"`).}
#'   \item{reference_2f}{`character`. Reference two-field allele for the P-group.}
#' }
#'
#' @details
#' P-group notation (e.g., `"A*01:01P"`) indicates that multiple alleles share
#' the same protein sequence in the antigen recognition domain. This table
#' maps P-groups to their representative reference alleles for serology lookup.
#'
#' @source WMDA nomenclature files from IMGT/HLA GitHub repository
#'   (\url{https://github.com/ANHIG/IMGTHLA/tree/Latest/wmda}).
#'
#' @usage data(deepMatchR_wmda_pgroups)
#' @seealso \code{\link{toSerology}}, \code{\link{deepMatchR_wmda_serology}}
#' @keywords datasets HLA WMDA P-group
"deepMatchR_wmda_pgroups"

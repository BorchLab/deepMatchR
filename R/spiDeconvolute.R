#' Deconvolute SAB and PRA findings
#'
#' Reconciles Single Antigen Bead (SAB) and Panel Reactive Antibody (PRA)
#' results for a single patient to identify concordant and discordant findings.
#'
#' @param SAB A data.frame of SAB results. Requires columns for antigen/allele
#'   and MFI (e.g., 'antigen', 'NormalValue').
#' @param PRA A data.frame of PRA results. Requires 'BeadID' and 'NormalValue'.
#' @param panelLong A long-format data.frame mapping antigens to beads in the
#'   PRA panel. Requires 'BeadID' and a key column (e.g., 'antigen').
#' @param epletMap (Optional) A data.frame for mapping antigens to eplets.
#'   Currently not implemented.
#' @param args A list of optional parameters to override defaults. See Details.
#'
#' @details
#' The `args` list can be used to customize the function's behavior.
#' Default values are:
#' \itemize{
#'   \item `keyCol`: "antigen"
#'   \item `praMfiCutoff`: 1000
#'   \item `sabMfiCutoff`: NULL (if set, overrides class-specific cutoffs. Can be a numeric value or a method like "robust")
#'   \item `sabMfiCutoffClassI`: 1500
#'   \item `sabMfiCutoffClassII`: 2500
#'   \item `k`: 5 (for robust cutoff method)
#'   \item `minSupportBeads`: 1
#'   \item `classFilter`: "both" (can be "I" or "II")
#' }
#'
#' The function performs the following steps:
#' 1. Determines PRA bead positivity using `praMfiCutoff`.
#' 2. Estimates an SAB MFI cutoff if not provided, using methods like "robust".
#' 3. Classifies SAB antigens as reactive based on the cutoff.
#' 4. Calculates PRA support for each antigen based on the panel composition.
#' 5. Partitions results into three tables: `concordant`, `sabOnly`, and `praOnly`.
#'
#' @return A list containing:
#' \itemize{
#'   \item `concordant`: A data.frame of antigens that are reactive in SAB and
#'     supported by PRA.
#'   \item `sabOnly`: A data.frame of antigens that are reactive in SAB but not
#'     supported by PRA.
#'   \item `praOnly`: A data.frame of antigens that are not reactive in SAB but
#'     are supported by PRA.
#'   \item `supportByBead`: A data.frame with detailed PRA support metrics for
#'     each antigen.
#'   \item `summary`: A list with summary counts and run information.
#'   \item `argsUsed`: A list of the final parameters used in the analysis.
#' }
#'
#' @examples
#' \dontrun{
#' # Example with mock data
#' sab <- data.frame(
#'   antigen = c("A1", "A2", "B7", "B8"),
#'   NormalValue = c(5000, 500, 3000, 200)
#' )
#' pra <- data.frame(
#'   BeadID = 1:4,
#'   NormalValue = c(2000, 800, 1500, 300)
#' )
#' panel <- data.frame(
#'   BeadID = c(1, 1, 2, 2, 3, 3, 4, 4),
#'   antigen = c("A1", "A2", "A1", "B7", "B7", "B8", "A2", "B8")
#' )
#'
#' result <- spiDeconvolute(sab, pra, panel)
#' result$concordant
#' result$summary
#' }
#'
#' @export
#' @importFrom data.table as.data.table copy fcase rbindlist
#' @importFrom utils modifyList
#'
spiDeconvolute <- function(SAB, PRA, panelLong, epletMap = NULL, args = list()) {

  # 1. Set up arguments
  default_args <- list(
    keyCol = "antigen",
    praMfiCutoff = 1000,
    sabMfiCutoff = NULL,
    sabMfiCutoffClassI = 1500,
    sabMfiCutoffClassII = 2500,
    k = 5,
    minSupportBeads = 1,
    classFilter = "both"
  )
  argsUsed <- modifyList(default_args, args)

  # Use data.table for performance
  sab_dt <- as.data.table(copy(SAB))
  pra_dt <- as.data.table(copy(PRA))
  panel_dt <- as.data.table(copy(panelLong))
  keyCol <- argsUsed$keyCol

  # 2. Pre-checks and Normalization
  if (!"class" %in% names(sab_dt)) {
      sab_dt[, class := fcase(
          grepl("^[ABCW]", sab_dt[[keyCol]]), "I",
          grepl("^[DRDPDQ]", sab_dt[[keyCol]]), "II",
          default = "Unknown"
      )]
  }

  if (argsUsed$classFilter != "both") {
      sab_dt <- sab_dt[class == argsUsed$classFilter]
      antigens_to_keep <- unique(sab_dt[[keyCol]])
      panel_dt <- panel_dt[get(keyCol) %in% antigens_to_keep]
  }

  # 3. Derive PRA bead positivity
  pra_pos <- .derivePraPositive(pra_dt, argsUsed$praMfiCutoff)

  # 4. Determine SAB reactivity
  sab_calls_list <- lapply(split(sab_dt, by = "class"), function(class_sab) {
      current_class <- class_sab$class[1]

      if (is.numeric(argsUsed$sabMfiCutoff)) {
        cutoff <- argsUsed$sabMfiCutoff
      } else if (is.character(argsUsed$sabMfiCutoff) && argsUsed$sabMfiCutoff == "robust") {
        cutoff <- .estimateSabCutoff(class_sab$NormalValue, k = argsUsed$k)
      } else {
        if (current_class == "I") {
            cutoff <- argsUsed$sabMfiCutoffClassI
        } else if (current_class == "II") {
            cutoff <- argsUsed$sabMfiCutoffClassII
        } else {
            cutoff <- argsUsed$sabMfiCutoffClassI # Fallback
        }
      }

      class_sab[, sab_mfi := NormalValue]
      class_sab[, sab_cutoff := cutoff]
      class_sab[, sab_reactive := sab_mfi >= sab_cutoff]
      return(class_sab[, c(keyCol, "sab_mfi", "sab_cutoff", "sab_reactive"), with = FALSE])
  })
  sab_calls <- rbindlist(sab_calls_list)

  # 5. Calculate agreement/support from PRA panel
  support_by_bead <- .agreementSupport(pra_pos, panel_dt, keyCol = keyCol)

  # 6. Partition results
  partitions <- .partitionResults(sab_calls, support_by_bead, argsUsed$minSupportBeads, keyCol = keyCol)

  # 7. Format Output
  summary_info <- list(
      concordant_n = nrow(partitions$concordant),
      sabOnly_n = nrow(partitions$sabOnly),
      praOnly_n = nrow(partitions$praOnly),
      total_sab_antigens = nrow(sab_calls),
      total_panel_antigens = length(unique(panel_dt[[keyCol]])),
      pra_positive_beads = sum(pra_pos$positive)
  )

  output <- list(
      concordant = as.data.frame(partitions$concordant),
      sabOnly = as.data.frame(partitions$sabOnly),
      praOnly = as.data.frame(partitions$praOnly),
      supportByBead = as.data.frame(support_by_bead),
      summary = summary_info,
      argsUsed = argsUsed
  )

  return(output)
}


# Helper functions for spiDeconvolute

#' @importFrom data.table as.data.table dcast
.incidenceFromPanel <- function(panelLong, keyCol = "antigen") {
  dt <- as.data.table(panelLong)
  # Ensure one entry per bead/key combination
  dt <- unique(dt, by = c("BeadID", keyCol))
  # Create the incidence matrix using dcast
  # value.var is implicitly 1s, fun.aggregate=length gives 1 if present, 0 if not
  incidence_matrix <- dcast(dt, as.formula(paste("BeadID ~", keyCol)),
                            fun.aggregate = length, value.var = keyCol)
  return(incidence_matrix)
}


#' @importFrom data.table as.data.table
.derivePraPositive <- function(PRA, praMfiCutoff) {
  dt <- as.data.table(PRA)
  if ("positive" %in% names(dt)) {
    dt[, positive := as.logical(positive)]
  } else {
    dt[, positive := NormalValue >= praMfiCutoff]
  }
  return(dt)
}

#' @importFrom stats median mad
.estimateSabCutoff <- function(sabMfi, method = "robust", k = 5) {
  if (method == "robust") {
    log_mfi <- log1p(sabMfi)
    med <- median(log_mfi, na.rm = TRUE)
    m <- mad(log_mfi, na.rm = TRUE)
    cutoff_log <- med + k * m
    cutoff <- expm1(cutoff_log)
    return(cutoff)
  }
  # Placeholder for other methods
  stop("Method '", method, "' not yet implemented.")
}


#' @importFrom data.table as.data.table
.agreementSupport <- function(praBeads, panelLong, keyCol = "antigen") {
  pra_dt <- as.data.table(praBeads)
  panel_dt <- as.data.table(panelLong)

  # Get positive beads
  positive_beads <- pra_dt[positive == TRUE, .(BeadID)]

  # Count total beads per antigen from panel
  total_beads_per_antigen <- panel_dt[, .N, by = keyCol]
  setnames(total_beads_per_antigen, "N", "totalBeads")

  # Count positive beads per antigen
  # 1. Filter panel for antigens on positive beads
  panel_on_pos_beads <- panel_dt[BeadID %in% positive_beads$BeadID]
  # 2. Count occurrences of each antigen
  pos_beads_per_antigen <- panel_on_pos_beads[, .N, by = keyCol]
  setnames(pos_beads_per_antigen, "N", "positiveBeads")

  # Merge total and positive counts
  support <- merge(total_beads_per_antigen, pos_beads_per_antigen, by = keyCol, all.x = TRUE)
  support[is.na(positiveBeads), positiveBeads := 0]

  # Calculate support fraction
  support[, supportFraction := positiveBeads / totalBeads]

  # Get list of supporting beads
  supporting_beads_list <- panel_on_pos_beads[, .(supportingBeads = list(unique(BeadID))), by = keyCol]
  support <- merge(support, supporting_beads_list, by = keyCol, all.x = TRUE)

  return(support)
}

#' @importFrom data.table as.data.table
.partitionResults <- function(sabCalls, agreementScores, minSupportBeads, keyCol = "antigen") {
  sab_dt <- as.data.table(sabCalls)
  agreement_dt <- as.data.table(agreementScores)

  # Merge SAB calls with agreement scores
  merged_data <- merge(sab_dt, agreement_dt, by = keyCol, all = TRUE)

  # Apply decision rules
  merged_data[, concordant := sab_reactive == TRUE & positiveBeads >= minSupportBeads]
  merged_data[, sabOnly := sab_reactive == TRUE & (positiveBeads < minSupportBeads | is.na(positiveBeads))]
  merged_data[, praOnly := (sab_reactive == FALSE | is.na(sab_reactive)) & positiveBeads >= minSupportBeads]

  # Create tidy tables
  concordant_tbl <- merged_data[concordant == TRUE]
  sabOnly_tbl <- merged_data[sabOnly == TRUE]
  praOnly_tbl <- merged_data[praOnly == TRUE]

  return(list(
    concordant = concordant_tbl,
    sabOnly = sabOnly_tbl,
    praOnly = praOnly_tbl
  ))
}

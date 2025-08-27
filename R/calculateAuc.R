#' Calculate Antigen AUC Based on MFI
#'
#' @description
#' This function provides a unified interface for analyzing Single Antigen Bead
#' (SAB) data against different antigenic features, such as eplets, 
#' cross-reactive groups (CREGs), or serology. It computes the proportion of features
#' that are positive above a sequence of MFI cut-offs and integrates that to 
#' obtain an area-under-the-curve (AUC).
#'
#' Depending on user arguments, it can either generate a ggplot or return
#' the AUC results as a data.frame.
#'
#' @param result_file A data frame of SAB results or a path to a CSV/XLS/XLSX file.
#' @param analysis_type Character. The type of analysis to perform.
#'   Must be one of `"eplet"`, `"creg"`, `"serology"`. 
#' @param feature_filter Integer. The minimum number of beads/alleles that must
#'   carry the feature for it to be included in the analysis. Defaults to `3`.
#' @param percPos_filter Numeric (0-1). The minimum proportion of beads that must
#'   be positive for at least one cut-off to keep the feature. Defaults to `0.8`.
#' @param group_by Character. The aesthetic used to color curves in the plot.
#' @param label Logical. If `TRUE`, labels the curves on the plot. Defaults to `TRUE`.
#' @param cut_min,cut_max,cut_step Range and step for MFI cut-offs.
#' @param plot_results Logical. If `TRUE` (default), returns a `ggplot` object.
#'   If `FALSE`, returns a summarized tibble with AUC results.
#' @param palette Character. A color palette name (see `grDevices::hcl.pals`)
#'   or a custom palette function. Defaults to `"spectral"`.
#' @param evidence_level For eplet analysis, a character vector of evidence levels to keep.
#' @param eplet_filter For eplet analysis, the minimum number of times an eplet must appear.
#' @param top_eplets For eplet analysis, the maximum number of top eplets to display.
#' @param creg_filter For CREG analysis, the minimum number of times a CREG must appear.
#' @param serology_filter For serology analysis, a filter to be applied.
#' @param ... Additional arguments passed to the plot theme.
#'
#' @return Either a `ggplot` object or a data.frame with AUC results. The data.frame
#'   will contain columns for the feature (`eplet`, `creg`, `serology`), `AUC`,
#'   `norm_AUC`, `total_count`, and `loci`.
#'
#' @importFrom ggplot2 ggplot aes geom_line xlim ylim labs scale_color_manual
#' @importFrom directlabels geom_dl last.points
#' @importFrom pracma trapz
#' @export
calculateAUC <- function(result_file,
                         analysis_type,
                         feature_filter = 3,
                         percPos_filter = 0.8,
                         group_by = NULL,
                         label = TRUE,
                         cut_min = 250,
                         cut_max = 10000,
                         cut_step = 250,
                         plot_results = TRUE,
                         palette = "spectral",
                         evidence_level = c("A1", "A2", "B", "D"),
                         eplet_filter = 3,
                         top_eplets = 10,
                         creg_filter = 3,
                         serology_filter = NULL,
                         ...) {

  
  # --- 1. Configure analysis based on type ---
  if (tolower(analysis_type) == "eplet") {
    config <- list(
      feature_col = "eplet",
      data = data.table::as.data.table(deepMatchR::deepMatchR_eplets),
      evidence_level = evidence_level,
      top_eplets = top_eplets,
      default_group_by = "eplet"
    )
  } else if (tolower(analysis_type) == "creg") {
    config <- list(
      feature_col = "CREG",
      data = data.table::as.data.table(deepMatchR::deepMatchR_cregs),
      default_group_by = "CREG"
    )
  } else if (tolower(analysis_type) == "serology") {
    config <- list(
      feature_col = "serology",
      data = data.table::as.data.table(deepMatchR::deepMatchR_cregs),
      default_group_by = "serology"
    )
  } else {
    stop("`analysis_type` must be one of 'eplet', 'creg' or 'serology'.")
  }
  
  if (is.null(group_by)) {
    group_by <- config$default_group_by
  }
  
  # --- 2. Load and process SAB data ---
  if (is.character(result_file)) {
    result0 <- .loadData(result_file)
  } else {
    result0 <- result_file
  }
  .checkSAB(result0)
  result <- data.table::as.data.table(.processSAB(result0))
  
  # --- 3. Create combinations of alleles and MFI cutoffs ---
  cutoffs <- seq(cut_min, cut_max, cut_step)
  class_alleles <- unique(result[, .(allele, mfi_min)])
  
  summary_dt <- data.table::CJ(allele = unique(class_alleles$allele), cut = cutoffs)
  summary_dt <- merge(summary_dt, class_alleles, by = "allele", all.x = TRUE)
  summary_dt <- summary_dt[mfi_min > cut, .(allele, cut)]
  
  # --- 4. Prepare feature dictionary (Eplet or CREG) ---
  feature_data <- config$data[allele %in% class_alleles$allele]
  
  if (analysis_type == "eplet" && !is.null(config$evidence_level)) {
    feature_data <- feature_data[evidence %in% config$evidence_level]
    if(nrow(feature_data) == 0) {
      stop("`evidence_level` filtering criteria did not produce any results.")
    }
  }
  
  feature_data[, count := .N, by = c(config[["feature_col"]], "allele")]
  feature_data[, subtotal := .N, by = c(config[["feature_col"]])]
  
  # --- 5. Calculate proportion positive for each feature × cut-off pair ---
  analysis_dt <- merge(summary_dt, feature_data, by = "allele", allow.cartesian = TRUE)
  analysis_dt[, loci := sub("\\*.*", "", allele)]
  analysis_dt <- analysis_dt[!is.na(cut)]

  analysis_dt[, positive_count := sum(count, na.rm = TRUE), by = c(config[["feature_col"]], "cut")]
  analysis_dt[, percent_positive := positive_count / subtotal]
  analysis_dt[, pp_max := max(percent_positive, na.rm = TRUE), by = c(config[["feature_col"]])]
  data.table::setorder(analysis_dt, -subtotal, -percent_positive)
  
  # --- 6. Apply user filters ---
  if (!is.null(feature_filter))
    analysis_dt <- analysis_dt[subtotal >= feature_filter]
  
  if (!is.null(percPos_filter))
    analysis_dt <- analysis_dt[pp_max >= percPos_filter]
  
  analysis_dt[, loci := paste0(unique(loci), collapse = "; "), by = c(config[["feature_col"]])]
  
  # --- 7. Calculate AUC ---
  feature_AUC <- analysis_dt[, .(
    AUC = pracma::trapz(cut, percent_positive),
    total_count = unique(subtotal)[1],
    loci = paste0(unique(loci), collapse = "; ")
  ), by = c(config[["feature_col"]])]
  feature_AUC[, norm_AUC := AUC / cut_max]
  
  # --- 8. Generate Plot or Return data.table ---
  if (plot_results) {
    if(group_by == "evidence_level") {
      group_by <- "evidence"
    }
    
    plot_data <- analysis_dt
    
    if (analysis_type == "eplet" && !is.null(config$top_eplets)) {
      top_features_vec <- feature_AUC[order(-norm_AUC)][1:config$top_eplets, get(config$feature_col)]
      plot_data <- plot_data[get(config$feature_col) %in% top_features_vec]
    }
    
    p <- ggplot(plot_data, aes(x = cut, y = percent_positive,
                               colour = .data[[group_by]], 
                               group = .data[[config$feature_col]])) +
      geom_line() +
      xlim(0, ifelse(label, cut_max + 1500, cut_max)) +
      ylim(0, 1) +
      .themeMatchR(...) + 
      labs(x = "Cut-off (MFI)", y = "Proportion Positive") +
      scale_color_manual(
        values = .colorizer(palette, length(unique(plot_data[[group_by]])))
      )
    
    if (label) {
      p <- p + geom_dl(aes(label = .data[[config$feature_col]]),
                       method = list("last.points", cex = 0.8))
    }
    return(p)
    
  } else {
    return(feature_AUC)
  }
}


# --- ALIAS WRAPPER FUNCTIONS ---

#' @rdname calculateAUC
#' @export
epletAUC <- function(result_file,
                     evidence_level = c("A1", "A2"),
                     eplet_filter = 3,
                     top_eplets = 10,
                     ...) {
  calculateAUC(
    result_file = result_file,
    analysis_type = "eplet",
    feature_filter = eplet_filter,
    evidence_level = evidence_level,
    top_eplets = top_eplets,
    ...
  )
}

#' @rdname calculateAUC
#' @export
cregAUC <- function(result_file,
                    creg_filter = 3,
                    ...) {
  calculateAUC(
    result_file = result_file,
    analysis_type = "creg",
    feature_filter = creg_filter,
    ...
  )
}

#' @rdname calculateAUC
#' @export
serologyAUC <- function(result_file,
                        serology_filter = 3,
                         ...) {
  calculateAUC(
    result_file = result_file,
    analysis_type = "serology",
    feature_filter = serology_filter,
    ...
  )
}
#' Plot Eplet Results from SPI Assay
#'
#' This function reads in SAB (single antigen bead) results (either as a data frame or
#' as a file path to a CSV/XLS/XLSX file), processes the data to quantify eplet-specific
#' positivity based on a specified MFI cutoff, and then generates one of three plot types:
#' a treemap, a bar plot, or an AUC plot. The eplet annotations are joined from an internal
#' database, and the resulting plot is colored by evidence level.
#'
#' @param result_file A data frame containing SAB results or a character string specifying
#'   the path to a SAB file in CSV, XLS, or XLSX format.
#' @param cutoff Numeric. Threshold for MFI to separate positive and negative beads.
#'   Default is 2000.
#' @param evidence_level Character vector indicating the antibody reactivity levels to keep.
#'   Defaults to \code{c("A1", "A2")}, which represent antibody-confirmed eplets.
#'   Other acceptable levels include \code{"B"}, \code{"D"}, or \code{NULL} to apply no filter.
#' @param group_by A character string or indicating the coloring grouping for 
#'   the plot, default is `eplet`. Other options include `loci` or 
#'   `evidence_level`.
#' @param eplet_filter Integer. (Only used when \code{plot_type = "AUC"}) 
#'   Specifying the minimum number of times an eplet must appear in the assay 
#'   before calculating the AUC. Defaults to `3`.
#' @param percPos_filter Numeric value. (Only used when \code{plot_type = "AUC"}) 
#'   Value between 0 and 1 specifying the minimum relative proportion of 
#'   positive beads (per eplet) to include in the final summary. Eplets with a 
#'   proportion below this threshold are excluded. Defaults to 0.4.
#' @param plot_type Character. Type of plot to generate. Must be one of \code{"treemap"},
#'   \code{"bar"}, or \code{"AUC"}. Defaults to \code{"treemap"}.
#' @param cut_min Integer. (Only used when \code{plot_type = "AUC"}) The minimum MFI cutoff value.
#'   Defaults to 250.
#' @param cut_max Integer. (Only used when \code{plot_type = "AUC"}) The maximum MFI cutoff value.
#'   Defaults to 10000.
#' @param cut_step Integer. (Only used when \code{plot_type = "AUC"}) The increment step between
#'   the minimum and maximum MFI cutoff values. Defaults to 250.
#' @param top_eplets Integer. The maximum number of top eplets to display in the bar or AUC plot.
#'   Defaults to 10.
#' @param palette Character. A color palette name (from \link[grDevices]{hcl.pals}) or a custom
#'   palette function to use for the plot. Defaults to \code{"spectral"}.
#' @param ... Additional arguments passed to the ggplot theme
#'
#' @return A \code{ggplot} object visualizing eplet counts (or AUC values) according to the
#'   specified parameters.
#' @importFrom data.table setkey 
#' @examples
#' # Using a data frame:
#' plotEplets(deepMatchR_example[[1]], 
#'            cutoff = 2000, 
#'            evidence_level = c("A1", "A2", "B"),
#'            percPos_filter = 0.4, 
#'            plot_type = "treemap")
#' @export
plotEplets <- function(result_file,
                       cutoff = 2000,
                       evidence_level = c("A1", "A2"),
                       group_by = "eplet",
                       eplet_filter = 3,
                       percPos_filter = 0.4,
                       plot_type = c("treemap", "bar", "AUC"),
                       cut_min = 250,
                       cut_max = 10000,
                       cut_step = 250,
                       top_eplets = 10,
                       palette = "spectral",
                       ...) {
  if(group_by == "evidence_level") group_by <- "evidence"
  eplet_data <- data.table::as.data.table(deepMatchR::deepMatchR_eplets)
  plot_type <- match.arg(plot_type)
  
  result0 <- if (inherits(result_file, "character")) .loadData(result_file) else result_file
  .checkSAB(result0)
  
  if (plot_type == "AUC") {
    summary_dt <- epletAUC(result_file = result0,
                           evidence_level = evidence_level,
                           plot_results = FALSE,
                           eplet_filter = eplet_filter,
                           percPos_filter = percPos_filter,
                           cut_min = cut_min,
                           cut_max = cut_max,
                           cut_step = cut_step)
    
    setkey(eplet_data, eplet)
    summary_dt <- eplet_data[summary_dt, on = "eplet", mult = "first"]
  } else {
    result <- data.table::as.data.table(.processSAB(result0))
    
    ep_analysis <- merge(result, eplet_data, by = "allele", all.x = TRUE, allow.cartesian=TRUE)
    ep_analysis[, loci := sub("\\*.*", "", allele)]
    
    if (!is.null(evidence_level)) {
      ep_analysis <- ep_analysis[evidence %in% evidence_level]
    }
    
    ep_analysis[, positive.bead := ifelse(NormalValue >= cutoff, 1, 0)]

    summary_dt <- ep_analysis[, .(
      loci = paste0(unique(loci), collapse = "; "),
      count_above = sum(positive.bead),
      count_total = .N,
      evidence_level = unique(evidence)
    ), by = eplet]

    summary_dt[, pp_max := round(count_above / count_total, 2)]
    
    if (!is.null(percPos_filter)) {
      summary_dt <- summary_dt[pp_max >= percPos_filter]
    }
  }
  
  summary_dt[[group_by]] <- factor(summary_dt[[group_by]],
                                   levels = .alphanumericalSort(unique(summary_dt[[group_by]])))

  color.palette <- .colorizer(palette, length(unique(summary_dt[[group_by]])))
  
  if (plot_type == "treemap") {
    plot <- ggplot(summary_dt, aes(area = abs(count_above) * pp_max,
                                   fill = .data[[group_by]],
                                   label = eplet,
                                   subgroup = loci)) +
      geom_treemap() +
      geom_treemap_text(aes(label = paste(eplet, "\n", pp_max * 100, "%\n",
                                          count_above, "of", count_total)),
                        place = "centre", grow = FALSE, min.size = 1, color = "black") +
      geom_treemap_subgroup_border(color = "black", size = 2) +
      geom_treemap_subgroup_text(place = "centre", grow = TRUE, alpha = 0.3, colour = "black") +
      scale_fill_manual(values = color.palette) +
      .themeMatchR(grid_lines = "No") +
      labs(fill = group_by) +
      theme(plot.background = element_blank())
    
  } else if (plot_type == "bar") {
    y.label <- "Proportion of Positive Beads"
    y <- "pp_max"
    data.table::setorder(summary_dt, -pp_max)
    ranked_data <- summary_dt[1:min(top_eplets, .N)]
    
    plot <- ggplot(ranked_data, aes(x = reorder(eplet, pp_max), y = .data[[y]],
                                    fill = .data[[group_by]])) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.25) +
      coord_flip(clip = "off") +
      geom_text(aes(label = loci), size = 2, hjust = -0.05) +
      labs(fill = group_by, y = y.label) +
      .themeMatchR() +
      theme(axis.title.y = element_blank()) +
      scale_fill_manual(values = color.palette)
    
  } else { # plot_type == "AUC"
    y.label <- "Normalized AUC"
    y <- "norm_AUC"
    data.table::setorder(summary_dt, -norm_AUC)
    ranked_data <- summary_dt[1:min(top_eplets, .N)]
    
    label.max <- round(max(ranked_data$norm_AUC), 2) - 0.05
    
    plot <- ggplot(ranked_data, aes(x = reorder(eplet, norm_AUC), y = .data[[y]],
                                    fill = .data[[group_by]])) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.25) +
      coord_flip(clip = "off") +
      geom_text(aes(label = loci, hjust = ifelse(.data[[y]] > label.max, 1.1, -0.1)), size = 2) +
      labs(fill = group_by, y = y.label) +
      .themeMatchR() +
      theme(axis.title.y = element_blank()) +
      scale_fill_manual(values = color.palette)
  }
  
  return(plot)
}

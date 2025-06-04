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
#'
#' @return A \code{ggplot} object visualizing eplet counts (or AUC values) according to the
#'   specified parameters.
#' @examples
#' # Using a data frame:
#' plotEplets(deepMatchR_example[[1]], 
#'            cutoff = 2000, 
#'            evidence_level = c("A1", "A2", "B"),
#'            percPos_filter = 0.4, 
#'            plot_type = "treemap")
#' @importFrom stringr str_sort
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
                       palette = "spectral") {
  
  # Standardize plot type argument early
  plot_type <- match.arg(plot_type)
  
  # Load required eplet database
  data(deepMatchR_eplets)
  
  # 1. Read in data (data frame or file path)
  if (inherits(result_file, "character")) {
    result0 <- .loadData(result_file)
  } else {
    result0 <- result_file
  }
  
  # 2. Check for required SAB columns
  .checkSAB(result0)
  
  # 3. Process data based on the plot type
  if (plot_type == "AUC") {
    # For AUC plot, calculate AUC values using the internal function
    summary_df <- epletAUC(result_file = result0,
                           evidence_level = evidence_level,
                           plot_results = FALSE,
                           eplet_filter = eplet_filter,
                           percPos_filter = percPos_filter,
                           cut_min = cut_min,
                           cut_max = cut_max,
                           cut_step = cut_step)
    
    # Filter by evidence level if specified
    if (!is.null(evidence_level)) {
      summary_df <- summary_df[which(summary_df[["evidence_level"]] %in% evidence_level),]
    }
    
  } else {
    # Clean and process the SAB results
    result <- .processSAB(result0)
    
    # 4. Combine eplet annotations with assay data
    ep_analysis <- result %>%
      left_join(deepMatchR_eplets, 
                by = "allele", 
                relationship = "many-to-many") %>%
      mutate(loci = str_extract(allele, "^[^*]+"))
    
    # Filter by evidence level if specified
    if (!is.null(evidence_level)) {
      ep_analysis <- ep_analysis[which(ep_analysis[["evidence_level"]] %in% evidence_level),]
    }
    
    # 5. Quantify positive beads for each eplet and summarize
    summary_df <- ep_analysis %>%
      mutate(positive.bead = ifelse(NormalValue >= cutoff, 1, 0)) %>%
      group_by(eplet) %>%
      summarise(loci = str_c(unique(loci), collapse = "; "), 
                count_above = sum(positive.bead),
                count_total = n(), 
                pp_max = round(count_above / count_total, 2), 
                evidence_level = unique(evidence_level),
                .groups = "drop")
    
    # Filter summarized data based on percent positive
    if (!is.null(percPos_filter)) {
      summary_df <- summary_df %>%
        filter(pp_max >= percPos_filter)
    }
  }
  
  summary_df[[group_by]] <- factor(summary_df[[group_by]], 
                                   levels = str_sort(unique(summary_df[[group_by]]), numeric = TRUE))
  # Generate the color palette using internal helper function
  color.palette <- .colorizer(palette, length(unique(summary_df[[group_by]])))
  
  # 6. Create the requested plot type
  if (plot_type == "treemap") {
    plot <- ggplot(summary_df, aes(area = abs(count_above) * pp_max, 
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
      .dmrTheme() + 
      labs(fill = group_by) +
      theme(plot.background = element_blank())
    
  } else if (plot_type == "bar") {  
    # Prepare data for bar plot: show top eplets ranked by proportion positive
    y.label <- "Proportion of Positive Beads"
    y <- "pp_max"
    ranked_data <- summary_df %>%
      arrange(desc(pp_max)) %>%
      mutate(rank = row_number()) %>%
      filter(rank <= top_eplets)
    
    color.palette <- .colorizer(palette, length(unique(ranked_data[[group_by]])))
    
    plot <- ggplot(ranked_data, aes(x = reorder(eplet, desc(rank)), y = .data[[y]], 
                                    fill = .data[[group_by]])) +
      geom_bar(stat = "identity", color = "black", size = 0.25) +
      coord_flip(clip = "off") +
      geom_text(aes(label = loci), size = 2, hjust = -0.05) +
      labs(fill = group_by, y = y.label) +
      .dmrTheme() + 
      theme(axis.title.y = element_blank()) + 
      scale_fill_manual(values = color.palette)
    
  } else {  # plot_type == "AUC"
    # Prepare data for AUC plot: show top eplets ranked by normalized AUC
    y.label <- "Normalized AUC"
    y <- "norm_AUC"
    ranked_data <- summary_df %>%
      arrange(desc(norm_AUC)) %>%
      mutate(rank = row_number()) %>%
      filter(rank <= top_eplets)
    
    color.palette <- .colorizer(palette, length(unique(ranked_data[[group_by]])))
    
    label.max <- round(max(ranked_data$norm_AUC),2) - 0.05
    
    plot <- ggplot(ranked_data, aes(x = reorder(eplet, dplyr::desc(rank)), y = .data[[y]], 
                                    fill = .data[[group_by]])) +
      geom_bar(stat = "identity", color = "black", size = 0.25) +
      coord_flip(clip = "off") +
      geom_text(aes(label = loci, hjust = ifelse(.data[[y]] > label.max, 1.1, -0.1)), size = 2) + 
      labs(fill = group_by, y = y.label) +
      .dmrTheme(base_size = 10) + 
      theme(axis.title.y = element_blank()) + 
      scale_fill_manual(values = color.palette)
  }
  
  return(plot)
}

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
#' @param percPos_filter Numeric value between 0 and 1 specifying the minimum relative
#'   proportion of positive beads (per eplet) to include in the final summary.
#'   Eplets with a proportion below this threshold are excluded. Defaults to 0.4.
#' @param plot.type Character. Type of plot to generate. Must be one of \code{"treemap"},
#'   \code{"bar"}, or \code{"AUC"}. Defaults to \code{"treemap"}.
#' @param cut_min Integer. (Only used when \code{plot.type = "AUC"}) The minimum MFI cutoff value.
#'   Defaults to 250.
#' @param cut_max Integer. (Only used when \code{plot.type = "AUC"}) The maximum MFI cutoff value.
#'   Defaults to 10000.
#' @param cut_step Integer. (Only used when \code{plot.type = "AUC"}) The increment step between
#'   the minimum and maximum MFI cutoff values. Defaults to 250.
#' @param top.eplets Integer. The maximum number of top eplets to display in the bar or AUC plot.
#'   Defaults to 20.
#' @param palette Character. A color palette name (from \link[grDevices]{hcl.pals}) or a custom
#'   palette function to use for the plot. Defaults to \code{"inferno"}.
#'
#' @return A \code{ggplot} object visualizing eplet counts (or AUC values) according to the
#'   specified parameters.

#' @examples
#' # Using a data frame:
#' plotEplets(deepMatchR_example[[1]], 
#'            cutoff = 2000, 
#'            evidence_level = c("A1", "A2", "B"),
#'            percPos_filter = 0.4, 
#'            plot.type = "treemap")
#'
#' @export
plotEplets <- function(result_file,
                       cutoff = 2000,
                       evidence_level = c("A1", "A2"),
                       percPos_filter = 0.4,
                       plot.type = c("treemap", "bar", "AUC"),
                       cut_min = 250,
                       cut_max = 10000,
                       cut_step = 250,
                       top.eplets = 20,
                       palette = "spectral") {
  
  # Standardize plot type argument early
  plot.type <- match.arg(plot.type)
  
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
  if (plot.type == "AUC") {
    # For AUC plot, calculate AUC values using the internal function
    summary_df <- epletAUC(result_file = result0,
                           plot_results = FALSE,
                           percPos_filter = percPos_filter,
                           cut_min = cut_min,
                           cut_max = cut_max,
                           cut_step = cut_step)
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
      ep_analysis <- ep_analysis %>%
        filter(antibody_reactivity %in% evidence_level)
    }
    
    # 5. Quantify positive beads for each eplet and summarize
    summary_df <- ep_analysis %>%
      mutate(positive.bead = ifelse(NormalValue >= cutoff, 1, 0)) %>%
      group_by(epitope) %>%
      summarise(loci = str_c(unique(loci), collapse = "; "), 
                count_above = sum(positive.bead),
                count_total = n(), 
                pp_max = round(count_above / count_total, 2), 
                antibody_reactivity = unique(antibody_reactivity),
                .groups = "drop")
    
    # Filter summarized data based on percent positive
    if (!is.null(percPos_filter)) {
      summary_df <- summary_df %>%
        filter(pp_max >= percPos_filter)
    }
  }
  
  # Generate the color palette using internal helper function
  color.palette <- .colorizer(palette = palette, 
                              n = length(unique(summary_df$antibody_reactivity)))
  
  # 6. Create the requested plot type
  if (plot.type == "treemap") {
    plot <- ggplot(summary_df, aes(area = abs(count_above) * pp_max, 
                                   fill = antibody_reactivity, 
                                   label = epitope, 
                                   subgroup = loci)) +
      geom_treemap() +
      geom_treemap_text(aes(label = paste(epitope, "\n", pp_max * 100, "%\n", 
                                          count_above, "of", count_total)),
                        place = "centre", grow = FALSE, min.size = 1, color = "black") +
      geom_treemap_subgroup_border(color = "black", size = 2) +
      geom_treemap_subgroup_text(place = "centre", grow = TRUE, alpha = 0.3, colour = "black") +
      scale_fill_manual(values = color.palette) +
      theme_clean() +
      labs(title = "Eplet Counts", fill = "Evidence") +
      theme(plot.background = element_blank())
    
  } else if (plot.type == "bar") {  
    # Prepare data for bar plot: show top eplets ranked by proportion positive
    y.label <- "Proportion of Positive Beads"
    y <- "pp_max"
    ranked_data <- summary_df %>%
      arrange(desc(pp_max)) %>%
      mutate(rank = row_number()) %>%
      filter(rank <= top.eplets)
    
    plot <- ggplot(ranked_data, aes(x = reorder(epitope, desc(rank)), y = .data[[y]], 
                                    fill = antibody_reactivity)) +
      geom_bar(stat = "identity", color = "black", size = 0.25) +
      coord_flip(clip = "off") +
      geom_text(aes(label = loci), size = 2, hjust = -0.05) +
      labs(fill = "Evidence", y = y.label) +
      theme_clean() +
      theme(axis.title.y = element_blank()) + 
      scale_fill_manual(values = color.palette)
    
  } else {  # plot.type == "AUC"
    # Prepare data for AUC plot: show top eplets ranked by normalized AUC
    y.label <- "Normalized AUC"
    y <- "norm_AUC"
    ranked_data <- summary_df %>%
      arrange(desc(norm_AUC)) %>%
      mutate(rank = row_number()) %>%
      filter(rank <= top.eplets)
    
    plot <- ggplot(ranked_data, aes(x = reorder(epitope, desc(rank)), y = .data[[y]], 
                                    fill = antibody_reactivity)) +
      geom_bar(stat = "identity", color = "black", size = 0.25) +
      coord_flip(clip = "off") +
      geom_text(aes(label = loci), size = 2, hjust = -0.05) +
      labs(fill = "Evidence", y = y.label) +
      theme_clean() +
      theme(axis.title.y = element_blank()) + 
      scale_fill_manual(values = color.palette)
  }
  
  return(plot)
}

#' Plotting Eplet Results from SPI Assay
#'
#' @param result_file A data frame containing SAB results **or** a character 
#'   string specifying the path to a SAB file in CSV, XLS, or XLSX format.
#' @param cutoff Numeric. Threshold for MFI to separate postive and negative beads
#' @param evidence_level A character string or vector of character strings
#'   indicating the desired evidence levels to keep. Defaults to 
#'   `c("A1", "A2")`, representing antibody-confirmed eplets. 
#'   Other levels include: `B`, `D`, or `NULL` if no filter is desired.
#' @param percPos_filter A numeric value between 0 and 1 representing the 
#'   relative percent of positive beads (for a specific eplet) to use as a filter. 
#'   Eplets below this threshold are excluded. Defaults to `0.4`
#' @param plot.type Character. Type of plot to generate: "treemap", "bar" or "auc".
#' @param cut_min If plotting AUC, an integer specifying the minimum MFI cutoff 
#'   value. Defaults to `250`.
#' @param cut_max If plotting AUC, an integer specifying the maximum MFI cutoff 
#'   value. Defaults to `10000`.
#' @param cut_step If plotting AUC, an integer specifying the increment step 
#'   between the minimum and maximum MFI cutoff values. Defaults to `250`.
#' @param palette A character string indicating the color palette to use when 
#'   plotting. Should be one of the palettes available through 
#'   \link[grDevices]{hcl.pals} or a custom function. Defaults to `"inferno"`.
#' 
#' @return A ggplot object visualizing eplet counts.

plotEplets <- function(result_file,
                       cutoff = 2000,
                       evidence_level = c("A1", "A2"),
                       percPos_filter = 0.4,
                       cut_min = 250,
                       cut_max = 10000,
                       cut_step = 250,
                       top.eplets = 20,
                       palette = "inferno") {
  
  # Load required eplet database 
  data(deepMatchR_eplets)
  
  # 1. Read in data (data frame or file path)
  if (inherits(result_file, "character")) {
    result0 <- .loadData(result_file)
  } else {
    # If it's already a data frame
    result0 <- result_file
  }
  
  # 2. Check if the incoming data has the required SAB columns
  .checkSAB(result0)
  
  if(plot.type == "AUC") {
    summary_df <- epletAUC(result_file = result0,
                           plot_results = FALSE,
                           percPos_filter = percPos_filter,
                           cut_min = cut_min,
                           cut_max = cut_max,
                           cut_step = cut_step)
  } else {
    # 3. Clean up and organize the screening result
    result <- .processSAB(result0)
    
  
    # 4. Combining eplets with assay data
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
    
   # 5. Quantify positive beads
    summary_df <- ep_analysis %>%
       mutate(positive.bead = ifelse(NormalValue >= cutoff, 1, 0)) %>%
       group_by(epitope) %>%
       summarise(loci = str_c(unique(loci), collapse = "; "), 
                 count_above = sum(positive.bead),
                 count_total = n(), 
                 pp_max = round(count_above/count_total, 2), 
                 antibody_reactivity = unique(antibody_reactivity))
    
    #Filter based on percent positive
    if (!is.null(percPos_filter)) {
      ep_analysis <- ep_analysis %>%
        filter(pp_max >= percPos_filter)
    }
    
  }
  color.palette <- .colorizer(palette = palette, 
                              n=length(unique(summary_df$antibody_reactivity)))
  
  plot.type <- match.arg(plot.type)
  
  if (plot.type == "treemap") {
    plot <- ggplot(summary_df, aes(area = abs(count_above) * pp_max, 
                                     fill = antibody_reactivity, 
                                     label = epitope, 
                                     subgroup = loci)) +
      geom_treemap() +
      geom_treemap_text(aes(label = paste(epitope, "\n", pp_max * 100, "%\n", count_above, "of", count_total)),
                        place = "centre", grow = FALSE, min.size = 1, color = "black") +
      geom_treemap_subgroup_border(color = "black", size = 2) +
      geom_treemap_subgroup_text(place = "centre", grow = TRUE, alpha = 0.3, colour = "black") +
      scale_fill_manual(values = color.palette) +
      theme_clean() +
      labs(title = "Eplet Counts", fill = "Evidence") +
      theme(plot.background = element_blank())
    
  } else {
    if(plot.type == "bar") {  
      y.label = "Proportion of Positive Beads"
      y <- "pp_max"
      ranked_data <- summary_df %>%
        arrange(desc(pp_max)) %>%
        mutate(rank = row_number()) %>%
        filter(rank <= top.eplets)
      
    } else {
      y.label = "Normalized AUC"
      y <- "norm_AUC"
      ranked_data <- summary_df %>%
        arrange(desc(norm_AUC)) %>%
        mutate(rank = row_number()) %>%
        filter(rank <= top.eplets)
    }
   
    
    plot <- ggplot(ranked_data, aes(x = reorder(epitope, desc(rank)), y = .data[[y]], fill = antibody_reactivity)) +
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

  
  
  
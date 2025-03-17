#' Calculate Eplet AUC based on MFI
#'
#' @description
#' This function reads in the specified files for SAB data and a screening result,
#' cleans and organizes the data, calculates how many eplets are positive above
#' a range of MFI cutoffs, computes the percentage of positive eplets, filters
#' them based on user-defined criteria, and finally calculates the area under the
#' curve (AUC) for each eplet. Depending on user arguments, it can either
#' generate and print a ggplot or export the results as a CSV file.
#'
#' @param result_file A data frame with SAB results or a character string 
#' specifying the path SAB file in csv or excel format
#' @param evidence_level A character or vector of characters to asset the eplet 
#' evidence level. Default is antibody-confirmed `c("A1", "A2")`
#' @param eplet_filter An integer specifying the minimum number of times the eplet
#' must appear in the assay to calculate the AIC
#' @param percPos_filter A value between 0 and 1 that represent the relative 
#' percent of positive beads with the specific eplet to use as a filter. 
#' @param cut_min An integer specifying the minimum cutoff value. Defaults to 250.
#' @param cut_max An integer specifying the maximum cutoff value. Defaults to 10000.
#' @param cut_step An integer specifying the increment step between min and max 
#' cutoffs. Defaults to 250.
#' @param plot_results Logical. If `TRUE`, a ggplot object is generated and 
#' printed. Defaults to `TRUE`.
#' @param export_results Logical. If `TRUE`, a CSV file of results is exported. 
#' Defaults to `FALSE`.
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}.
#'
#' @return A tibble (data.frame) containing the eplet analysis, including the 
#' AUC column or a ggplot.
#' 
#' @examples
#' #Plot Data
#' epletAUC(deepMatchR_example[[1]],
#'          plot_results = TRUE, 
#'          percPos_filter = 0.9)
#'                    
#' #Export Data
#' result <- epletAUC(deepMatchR_example[[1]],
#'                    plot_results = FALSE,
#'                    export_results = TRUE, 
#'                    percPos_filter = 0.9)
#' @importFrom dplyr filter mutate select distinct arrange group_by ungroup 
#' summarise relocate bind_rows
#' @importFrom tidyr unnest_longer separate_longer_delim
#' @importFrom stringr str_detect str_extract str_replace_all
#' @importFrom ggplot2 ggplot aes geom_line geom_dl scale_colour_discrete xlim 
#' ylim theme labs
#' @importFrom directlabels dl.combine last.points
#' @importFrom janitor clean_names
#' @importFrom purrr map
#' @importFrom pracma trapz
#'
#' @export
epletAUC <- function(result_file,
                     evidence_level = c("A1", "A2"),
                     eplet_filter = 3,
                     percPos_filter = 0.8,
                     cut_min = 250,
                     cut_max = 10000,
                     cut_step = 250,
                     plot_results = TRUE,
                     export_results = FALSE, 
                     palette = "inferno") {
  
  data(deepMatchR_eplets)
  
  ## 1. Read in data files
  if(inherits(result_file, "character")) {
    result0   <- .loadData(result_file)
  } else {
    result0 <- result_file
  }
  # Check if required columns are present
  .checkSAB(result0)
  
  # 2. Clean up and organize the screening result
  result <- result0 %>%
    select(BeadID, SpecAbbr, Specificity, NormalValue) %>%
    distinct(Specificity, .keep_all = TRUE) %>%
    mutate(antigen = str_extract(SpecAbbr, '[ABCDRQP][:alnum:]+')) %>%
    mutate(bw46 = str_extract(SpecAbbr, 'Bw[46]')) %>%
    mutate(Specificity_truncated = str_extract(Specificity, '[ABCD].*[0-9]')) %>%
    mutate(allele = str_replace_all(Specificity_truncated, ",-,", "_")) %>%
    select(-SpecAbbr, -Specificity, -Specificity_truncated) %>%
    relocate(BeadID, antigen, bw46, allele, NormalValue) %>%
    separate_longer_delim(allele, "_") %>%  # split dimer beads into separate rows
    mutate(mfi_min = min(NormalValue), .by = allele) %>%  
    arrange(allele, desc(NormalValue)) %>%
    filter(!is.na(allele)) %>%
    distinct(allele, .keep_all = TRUE)     # keep top row for each allele
  
  # 3. Create all combinations of class alleles and cutoffs in a vectorized way
  cutoffs <- seq(cut_min, cut_max, cut_step)
  class_alleles <- result %>%
    select(allele, mfi_min)
  
  # Expand grid: each allele with each cutoff
  summary_df <- expand.grid(allele = class_alleles$allele, cut = cutoffs) %>%
    as_tibble() %>%
    left_join(class_alleles, by = "allele") %>%
    filter(mfi_min > cut) %>%
    select(allele, cut)
  
  # 4. Subset Eplet Dictionary for SAB-specific Information
  assay_alleles <- deepMatchR_eplets[deepMatchR_eplets$allele %in% class_alleles$allele,]
  
  if(!is.null(evidence_level)) {
    assay_alleles <- assay_alleles %>%
      filter(antibody_reactivity %in% evidence_level)
  }
  
  assay_alleles <- assay_alleles %>%
                    group_by(epitope, allele) %>%
                    mutate(count = n()) %>%  
                    ungroup() %>%
                    group_by(epitope) %>%
                    mutate(subtotal = n())
  
  # 5. Join with sab1_ep and analyze eplet reactivity
  ep_analysis <- summary_df %>%
                    left_join(assay_alleles, 
                              by = "allele", 
                              relationship = "many-to-many") %>%
    filter(!is.na(cut)) %>%
    group_by(epitope, cut) %>%
    mutate(
      positive_count   = sum(count),
      percent_positive = sum(count) / subtotal
    ) %>%
    group_by(epitope) %>%
    mutate(pp_max = max(percent_positive)) %>%
    arrange(desc(subtotal), desc(percent_positive)) %>%
    ungroup()
  
  if(!is.null(eplet_filter)) {
    ep_analysis <- ep_analysis %>%
                      filter(subtotal >= 3)
  }
  
  if(!is.null(percPos_filter)) {
    ep_analysis <- ep_analysis %>%
      filter(pp_max >= percPos_filter)
  }
  
  
  # 6. Generate a plot if requested
  if (plot_results) {
    # We can color by eplet; if the number of eplets is large, you may want
    # to choose a different approach or facet by eplet.
    plot <- ggplot(ep_analysis, aes(x = cut, y = percent_positive)) +
      geom_line(aes(color = epitope)) +
      scale_colour_discrete(guide = 'none') +
      xlim(0, 12000) +
      ylim(0, 1) +
      theme(aspect.ratio = 1) +
      geom_dl(aes(label = epitope), method = list(dl.combine("last.points")), cex = 0.8) +
      labs(x = "Cutoff (MFI)",
           y = "Proportion Positive") + 
      guides(color = "none") + 
      scale_color_manual(values = .colorizer(palette, length(unique(ep_analysis$epitope))))+ 
      theme_minimal()
    
    return(plot)
  } else {
  # 7. Return Data Frame Summary
    
    ep_auc <- ep_analysis %>%
      group_by(epitope) %>%
      summarize(
        AUC = trapz(x = cut, y = percent_positive),
        subtotal = first(subtotal)  
      ) %>%
      ungroup()
    
    return(ep_auc)
  }
  
}

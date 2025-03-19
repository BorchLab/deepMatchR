#' Calculate Eplet AUC Based on MFI
#'
#' @description
#' This function reads in Single Antigen Bead (SAB) data from either a data frame
#' or a file path (CSV, XLS, or XLSX). It cleans and organizes the data, 
#' calculates how many eplets are positive above a range of MFI cutoffs,
#' computes the percentage of positive eplets, optionally filters them based on 
#' user-defined criteria, and finally calculates the area under the curve (AUC)
#' for each eplet. Depending on user arguments, it can either generate and print 
#' a ggplot **or** return/export the AUC-based results as a CSV.
#'
#' @param result_file A data frame containing SAB results **or** a character 
#'   string specifying the path to a SAB file in CSV, XLS, or XLSX format.
#' @param evidence_level A character string or vector of character strings
#'   indicating the desired evidence levels to keep. Defaults to 
#'   `c("A1", "A2")`, representing antibody-confirmed eplets. 
#'   Other levels include: `B`, `D`, or `NULL` if no filter is desired.
#' @param eplet_filter An integer specifying the minimum number of times an 
#'   eplet must appear in the assay before calculating the AUC. Defaults to `3`.
#' @param percPos_filter A numeric value between 0 and 1 representing the 
#'   relative percent of positive beads (for a specific eplet) to use as a filter. 
#'   Eplets below this threshold are excluded. Defaults to `0.8`.
#' @param cut_min An integer specifying the minimum MFI cutoff value. 
#'   Defaults to `250`.
#' @param cut_max An integer specifying the maximum MFI cutoff value. 
#'   Defaults to `10000`.
#' @param cut_step An integer specifying the increment step between the minimum 
#'   and maximum MFI cutoff values. Defaults to `250`.
#' @param plot_results Logical. If `TRUE`, the function returns and prints a 
#'   \code{ggplot} object illustrating the proportion of positive eplets at 
#'   each cutoff. If `FALSE`, the function returns a summarized tibble.
#'   Defaults to `TRUE`.
#' @param palette A character string indicating the color palette to use when 
#'   plotting. Should be one of the palettes available through 
#'   \link[grDevices]{hcl.pals} or a custom function. Defaults to `"inferno"`.
#'
#' @return If \code{plot_results = TRUE}, a \code{ggplot} object is returned. 
#'   If \code{plot_results = FALSE}, a tibble (data frame) is returned. If 
#'   \code{export_results = TRUE} and \code{plot_results = FALSE}, the results 
#'   are also written to a CSV file named \code{"epletAUC_results.csv"} in 
#'   the working directory.
#' 
#' @examples
#' # Example 1: Plot data
#' epletAUC(
#'   result_file = deepMatchR_example[[1]],
#'   plot_results = TRUE, 
#'   percPos_filter = 0.9
#' )
#'
#' # Example 2: Export data
#' result <- epletAUC(
#'   result_file = deepMatchR_example[[1]],
#'   plot_results = FALSE,
#'   export_results = TRUE, 
#'   percPos_filter = 0.9
#' )
#'
#' @importFrom dplyr filter mutate select distinct arrange group_by ungroup 
#'   summarise relocate left_join n
#' @importFrom tidyr unnest_longer separate_longer_delim
#' @importFrom stringr str_detect str_extract str_replace_all
#' @importFrom ggplot2 ggplot aes geom_line scale_colour_discrete xlim ylim theme 
#'   labs theme_minimal guides scale_color_manual
#' @importFrom directlabels geom_dl dl.combine last.points
#' @importFrom janitor clean_names
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
                     palette = "spectral") {
  
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
  
  # 3. Clean up and organize the screening result
  result <- .processSAB(result0)
  
  # 4. Create all combinations of alleles and user-specified MFI cutoffs
  cutoffs <- seq(cut_min, cut_max, cut_step)
  class_alleles <- result %>%
    select(allele, mfi_min)
  
  # Generate a grid of each allele x cutoff
  summary_df <- expand.grid(allele = class_alleles$allele, cut = cutoffs) %>%
    as_tibble() %>%
    left_join(class_alleles, by = "allele") %>%
    filter(mfi_min > cut) %>%
    select(allele, cut)
  
  # 5. Subset Eplet Dictionary to only those alleles found in the SAB data
  assay_alleles <- deepMatchR_eplets[deepMatchR_eplets$allele %in% class_alleles$allele, ]
  
  # Filter by evidence level if specified
  if (!is.null(evidence_level)) {
    assay_alleles <- assay_alleles %>%
      filter(antibody_reactivity %in% evidence_level)
  }
  
  # For each epitope-allele pair, note how many times it appears
  assay_alleles <- assay_alleles %>%
    group_by(epitope, allele) %>%
    mutate(count = n()) %>%
    ungroup() %>%
    group_by(epitope) %>%
    mutate(subtotal = n()) %>%
    ungroup()
  
  # 6. Join the sab data with eplet dictionary to analyze eplet reactivity
  ep_analysis <- summary_df %>%
    left_join(assay_alleles, by = "allele", relationship = "many-to-many") %>%
    mutate(loci = str_extract(allele, "^[^*]+")) %>% 
    filter(!is.na(cut)) %>%
    group_by(epitope, cut) %>%
    mutate(
      positive_count   = sum(count, na.rm = TRUE),
      percent_positive = positive_count / subtotal
    ) %>%
    group_by(epitope) %>%
    mutate(pp_max = max(percent_positive, na.rm = TRUE)) %>%
    arrange(desc(subtotal), desc(percent_positive)) %>%
    ungroup()
  
  # 7. Apply user-specified filters:
  if (!is.null(eplet_filter)) {
    ep_analysis <- ep_analysis %>%
      filter(subtotal >= eplet_filter)
  }
  
  if (!is.null(percPos_filter)) {
    ep_analysis <- ep_analysis %>%
      filter(pp_max >= percPos_filter)
  }
  
  # 8. If the user wants to plot results, generate a ggplot
  if (plot_results) {
    
    plot <- ggplot(ep_analysis, aes(x = cut, y = percent_positive)) +
      geom_line(aes(color = epitope)) +
      scale_colour_discrete(guide = 'none') +
      xlim(0, 12000) +
      ylim(0, 1) +
      theme_minimal() +
      theme(aspect.ratio = 1) +
      geom_dl(aes(label = epitope), method = list(dl.combine("last.points"))) +
      labs(
        x = "Cutoff (MFI)",
        y = "Proportion Positive"
      ) + 
      guides(color = "none") + 
      scale_color_manual(
        values = .colorizer(palette, length(unique(ep_analysis$epitope)))
      )
    
    return(plot)
    
  } else {
    # 9. Otherwise, compute area under the curve (AUC) and return a tibble
    ep_auc <- ep_analysis %>%
      group_by(epitope) %>%
      summarize(
        AUC      = trapz(x = cut, y = percent_positive),
        norm_AUC = AUC/cut_max,
        total_count = first(subtotal),
        antibody_reactivity = unique(antibody_reactivity), 
        loci = str_c(unique(loci), collapse = "; ")
      ) %>%
      ungroup()
    
    return(ep_auc)
  }
}

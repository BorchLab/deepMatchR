epletAUC <- function(result_file,
                     evidence_level = c("A1", "A2"),
                     group_by = "eplet",
                     label = TRUE,
                     eplet_filter = 3,
                     percPos_filter = 0.8,
                     cut_min = 250,
                     cut_max = 10000,
                     cut_step = 250,
                     plot_results = TRUE,
                     top_eplets = 10,
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
    dplyr::group_by(allele) %>%
    summarise(mfi_min = median(mfi_min))
  #TODO Future think about allow for other ways to summarise alleles
  
  # Generate a grid of each allele x cutoff
  summary_df <- expand.grid(allele = class_alleles$allele, cut = cutoffs) %>%
    as_tibble() %>%
    left_join(class_alleles, by = "allele") %>%
    dplyr::filter(mfi_min > cut) %>%
    dplyr::select(allele, cut)
  
  # 5. Subset Eplet Dictionary to only those alleles found in the SAB data
  assay_alleles <- deepMatchR_eplets[deepMatchR_eplets$allele %in% class_alleles$allele, ]
  
  # Filter by evidence level if specified
  if (!is.null(evidence_level)) {
    assay_alleles <- assay_alleles[assay_alleles[["evidence_level"]] %in% evidence_level,]
  }
  
  # For each eplet-allele pair, note how many times it appears
  assay_alleles <- assay_alleles %>%
    group_by(eplet, allele) %>%
    mutate(count = n()) %>%
    ungroup() %>%
    group_by(eplet) %>%
    mutate(subtotal = n()) %>%
    ungroup()
  
  # 6. Join the sab data with eplet dictionary to analyze eplet reactivity
  ep_analysis <- summary_df %>%
    left_join(assay_alleles, by = "allele", relationship = "many-to-many") %>%
    mutate(loci = sub("\\*.*", "", allele)) %>% 
    dplyr::filter(!is.na(cut)) %>%
    group_by(eplet, cut) %>%
    mutate(
      positive_count   = sum(count, na.rm = TRUE),
      percent_positive = positive_count / subtotal
    ) %>%
    group_by(eplet) %>%
    mutate(pp_max = max(percent_positive, na.rm = TRUE)) %>%
    arrange(desc(subtotal), desc(percent_positive)) %>%
    ungroup()
  
  # 7. Apply user-specified filters:
  if (!is.null(eplet_filter)) {
    ep_analysis <- ep_analysis %>%
      dplyr::filter(subtotal >= eplet_filter)
  }
  
  if (!is.null(percPos_filter)) {
    ep_analysis <- ep_analysis %>%
      dplyr::filter(pp_max >= percPos_filter)
  }
  
  # Relevel Eplet Loci after all filter and calculations
  ep_analysis <- ep_analysis %>%
      group_by(eplet) %>%
      mutate(loci = paste(unique(loci), collapse = "; "))
  
  # 8 compute area under the curve (AUC) and return a tibble
  ep_AUC <- ep_analysis %>%
    group_by(eplet) %>%
    summarize(
      AUC      = trapz(x = cut, y = percent_positive),
      norm_AUC = AUC/cut_max,
      total_count = unique(subtotal)[1],
      evidence_level = unique(evidence_level), 
      loci = paste(unique(loci), collapse = "; ")
    ) %>%
    ungroup()
  
  #9. If the user wants to plot results, generate a ggplot
  if (plot_results) {  
    
    top_eplet_vec <- ep_AUC %>%
      slice_max(order_by = norm_AUC, n = top_eplets) %>%
      pull(eplet)
    
    plot <- ep_analysis %>%
      subset(eplet %in% top_eplet_vec) %>%
    ggplot(aes(x = cut, y = percent_positive)) +
      geom_line(aes(color = .data[[group_by]], group = eplet)) +
      xlim(0, ifelse(label, cut_max + 1500, cut_max)) +
      ylim(0, 1) +
      .themeMatchR() + 
      labs(
        x = "Cutoff (MFI)",
        y = "Proportion Positive"
      ) + 
      scale_color_manual(
        values = .colorizer(palette, length(unique(ep_analysis[[group_by]])))) + 
      if (label) list(geom_dl(aes(label = eplet), method = list("last.points", cex = 0.8))) else list()
    
    return(plot)
    
  } else { # 10. Otherwise, return AUC
    
    return(ep_AUC)
  }
}

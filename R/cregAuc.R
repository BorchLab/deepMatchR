cregAUC <- function(result_file,
                    group_by       = "creg",
                    label          = TRUE,
                    creg_filter    = 3,
                    percPos_filter = 0.8,
                    cut_min        = 250,
                    cut_max        = 10000,
                    cut_step       = 250,
                    plot_results   = TRUE,
                    palette        = "spectral") {
  
  # Load CREG dictionary delivered with deepMatchR
  data(deepMatchR_cregs)
  
  # 1. Read in data (data frame or file path)
  if (is.character(result_file)) {
    result0 <- .loadData(result_file)
  } else {
    result0 <- result_file
  }
  
  # 2. Check if the incoming data has the required SAB columns
  .checkSAB(result0)           
  result <- .processSAB(result0)  
  
  # 3. Create all combinations of alleles and user-specified MFI cutoffs
  cutoffs <- seq(cut_min, cut_max, cut_step)
  class_alleles <- result %>% dplyr::select(allele, mfi_min)
  
  summary_df <- expand.grid(allele = class_alleles$allele, cut = cutoffs) |>
    as_tibble() |>
    left_join(class_alleles, by = "allele") |>
    dplyr::filter(mfi_min > cut) |>
    dplyr::select(allele, cut)
  
  # 4. Subset CREG Dictionary to only those alleles found in the SAB data
  assay_creg <- deepMatchR_cregs[deepMatchR_cregs$allele %in%
                                   class_alleles$allele, ]
  
  # per-CREG bookkeeping
  assay_creg <- assay_creg |>
    group_by(creg, allele) |> mutate(count = n()) |> ungroup() |>
    group_by(creg)         |> mutate(subtotal = n()) |> ungroup()
  
  # 5.  Calculate proportion positive for each CREG × cut-off pair
  cr_analysis <- summary_df |>
    left_join(assay_creg, by = "allele", relationship = "many-to-many") |>
    mutate(loci = str_extract(allele, "^[^*]+")) |>
    dplyr::filter(!is.na(cut)) |>
    group_by(creg, cut) |>
    mutate(
      positive_count   = sum(count, na.rm = TRUE),
      percent_positive = positive_count / subtotal
    ) |>
    group_by(creg) |>
    mutate(pp_max = max(percent_positive, na.rm = TRUE)) |>
    arrange(desc(subtotal), desc(percent_positive)) |>
    ungroup()

  # 6.  User filters
  if (!is.null(creg_filter))
    cr_analysis <- cr_analysis |> dplyr::filter(subtotal >= creg_filter)
  
  if (!is.null(percPos_filter))
    cr_analysis <- cr_analysis |> dplyr::filter(pp_max >= percPos_filter)
  
  ## Collapse loci for labelling
  cr_analysis <- cr_analysis |>
    group_by(creg) |>
    mutate(loci = str_c(unique(loci), collapse = "; "))
  
  # 7. If the user wants to plot results, generate a ggplot
  if (plot_results) {
    p <- ggplot(cr_analysis,
                aes(x = cut, y = percent_positive,
                    colour = .data[[group_by]], group = creg)) +
      geom_line() +
      xlim(0, ifelse(label, cut_max + 1500, cut_max)) +
      ylim(0, 1) +
      .themeMatchR() +
      labs(x = "Cut-off (MFI)",
           y = "Proportion Positive") +
      scale_color_manual(
        values = .colorizer(palette,
                            length(unique(cr_analysis[[group_by]]))))
    
    if (label)
      p <- p + geom_dl(aes(label = creg),
                       method = list("last.points", cex = 0.8))
    
    return(p)
  }
  
  # Otherwise, compute area under the curve (AUC) and return a tibble
  cr_AUC <- cr_analysis |>
    group_by(creg) |>
    summarise(
      AUC         = trapz(cut, percent_positive),
      norm_AUC    = AUC / cut_max,
      total_count = unique(subtotal)[1],
      loci        = str_c(unique(loci), collapse = "; ")
    ) |>
    ungroup()
  
  cr_AUC
}

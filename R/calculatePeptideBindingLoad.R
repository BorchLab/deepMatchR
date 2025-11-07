#' Calculate Peptide Binding Load from Mismatched Positions
#'
#' @description
#' Identifies ALL donor alleles not present in recipient, generates peptides around
#' mismatch positions, and predicts MHC binding to ALL recipient alleles regardless 
#' of locus. Adaptively tests Class I peptides (8-11mers) against Class I alleles 
#' (A, B, C) and Class II peptides (12-25mers) against Class II alleles (DRB1, DQB1, etc.).
#' 
#' @param recipient_geno Recipient HLA genotype object
#' @param donor_geno Donor HLA genotype object
#' @param loci Character vector of loci to analyze (NULL for all)
#' @param ic50_threshold IC50 threshold in nM for binding (default 500)
#' @param return What to return: "summary", "detailed", "by_recipient_allele"
#' @param parallel Use parallel processing
#' @param n_cores Number of cores for parallel processing
#' @param hla_env Basilisk environment for MHCnuggets
#' @param cache_dir Optional directory for sequence caching
#' @param compare_to How to select recipient sequences for comparison:
#'   "all" = compare to all recipient alleles at same locus (default)
#'   "closest" = compare to most similar recipient allele at same locus
#'
#' @return Results dataframe with binding predictions
#'
#' @importFrom dplyr group_by summarize arrange
#' @export
calculatePeptideBindingLoad <- function(recipient_geno,
                                        donor_geno,
                                        loci = NULL,
                                        ic50_threshold = 500L,
                                        return = c("summary", "detailed", "by_recipient_allele"),
                                        parallel = TRUE,
                                        n_cores = 2,
                                        hla_env = deepmatchrEnv(),
                                        cache_dir = NULL,
                                        compare_to = c("all", "closest")) {
  
  return <- match.arg(return)
  compare_to <- match.arg(compare_to)
  
  # Validate genotypes
  validateHlaGeno(recipient_geno)
  validateHlaGeno(donor_geno)
  
  # Determine loci to analyze
  shared_loci <- intersect(recipient_geno$locus_present, donor_geno$locus_present)
  # Removing DRA Comparison
  shared_loci <- shared_loci[shared_loci != "DRA"]
  if (!is.null(loci)) shared_loci <- intersect(shared_loci, loci)
  if (length(shared_loci) == 0) stop("No shared loci between donor and recipient.")
  
  # Helper function to determine MHC class from allele name
  getAlleleClass <- function(allele) {
    # Class II: DRB1, DRB3, DRB4, DRB5, DQA1, DQB1, DPA1, DPB1
    if (grepl("^D[RPQ][AB][0-9]", allele)) {
      return("II")
    } else {
      # Class I: A, B, C, E, F, G
      return("I")
    }
  }
  
  #Helper for selecting Class II models
  get_classII_pairs <- function(alleles, pair_by = FALSE, sep = "-") {
    # --- Basic validation ---
    if (!is.character(alleles))
      stop("'alleles' must be a character vector.")
    if (is.null(names(alleles)) || any(is.na(names(alleles))))
      stop("'alleles' must be a *named* character vector (no NA names).")
    
    # --- Extract per-locus unique alleles (unnamed) ---
    grab <- function(locus) unique(unname(alleles[names(alleles) == locus]))
    dpa1 <- grab("DPA1")
    dpb1 <- grab("DPB1")
    dqa1 <- grab("DQA1")
    dqb1 <- grab("DQB1")
    drb1 <- grab("DRB1")
    
    # --- Helper to build pairs ---
    make_pairs <- function(a, b) {
      if (length(a) == 0 || length(b) == 0) return(character(0))
      if (pair_by) {
        n <- min(length(a), length(b))
        if (n == 0) return(character(0))
        unique(paste0(a[seq_len(n)], sep, b[seq_len(n)]))
      } else {
        unique(paste0(rep(a, each = length(b)), sep, b))
      }
    }
    
    # --- Build outputs with presence checks ---
    out <- character(0)
    
    if (length(dpa1) == 0 && length(dpb1) > 0)
      warning("No DPA1 alleles found; skipping DP pairing.")
    if (length(dpb1) == 0 && length(dpa1) > 0)
      warning("No DPB1 alleles found; skipping DP pairing.")
    if (length(dqa1) == 0 && length(dqb1) > 0)
      warning("No DQA1 alleles found; skipping DQ pairing.")
    if (length(dqb1) == 0 && length(dqa1) > 0)
      warning("No DQB1 alleles found; skipping DQ pairing.")
    if (length(drb1) == 0)
      warning("No DRB1 alleles found; none will be returned for DRB1.")
    
    dp_pairs <- make_pairs(dpa1, dpb1)
    dq_pairs <- make_pairs(dqa1, dqb1)
    
    # Add DRB1 as-is (individual alleles)
    out <- c(out, dp_pairs, dq_pairs, drb1)
    
    unique(out)
  }
  
  # Extract alleles
  get_alleles_for_loci <- function(geno, loci_vec) {
    pat <- paste0("^(", paste(loci_vec, collapse = "|"), ")_")
    idx <- grep(pat, colnames(geno$data))
    alleles <- unlist(geno$data[1, idx, drop = FALSE])
    alleles <- alleles[!is.na(alleles) & nzchar(alleles)]
    names(alleles) <- sub("_\\d+$", "", names(alleles))  # Keep locus info
    alleles
  }
  
  r_alleles <- get_alleles_for_loci(recipient_geno, shared_loci)
  d_alleles <- get_alleles_for_loci(donor_geno, shared_loci)
  
  # Step 1: Identify mismatched alleles (donor alleles NOT in recipient)
  mismatched_alleles <- setdiff(d_alleles, r_alleles)
  
  if (length(mismatched_alleles) == 0) {
    message("No mismatched alleles found between donor and recipient.")
    return(data.frame(
      total_mismatched_peptides = 0,
      binding_peptides = 0,
      binding_percentage = 0,
      ic50_threshold = ic50_threshold
    ))
  }
  
  message(sprintf("Found %d mismatched donor alleles: %s", 
                  length(mismatched_alleles), 
                  paste(mismatched_alleles, collapse = ", ")))
  message(sprintf("Testing binding to %d recipient alleles: %s", 
                  length(r_alleles), 
                  paste(r_alleles, collapse = ", ")))
  
  # Classify recipient alleles by MHC class
  recipient_class <- sapply(r_alleles, getAlleleClass)
  class_I_recipients <- r_alleles[recipient_class == "I"]
  class_II_recipients <- r_alleles[recipient_class == "II"]
  
  message(sprintf("  - %d Class I recipient alleles: %s",
                  length(class_I_recipients),
                  paste(class_I_recipients, collapse = ", ")))
  message(sprintf("  - %d Class II recipient alleles: %s",
                  length(class_II_recipients),
                  paste(class_II_recipients, collapse = ", ")))
  
  # Step 2: Get protein sequences for all alleles
  all_alleles <- unique(c(mismatched_alleles, r_alleles))
  prot_seq_map <- batchGetSequences(all_alleles, type = "PROT", 
                                    n_cores = n_cores, cache_dir = cache_dir)
  
  # Step 3: Generate peptides from mismatched alleles
  # Generate BOTH Class I (8-11) and Class II (12-25) peptides
  generatePeptidesFromMismatches <- function(donor_seq, recipient_seq, 
                                             donor_allele, recipient_allele,
                                             class_I_lengths = 8:11,
                                             class_II_lengths = 12:25) {
    # Use quantifyMismatch to align and find mismatch positions
    mismatch_detail <- quantifyMismatch(
      donor_seq, 
      recipient_seq, 
      return = "detail"
    )
    
    # Get positions where there are mismatches
    mismatch_positions <- mismatch_detail$alignment_position[mismatch_detail$is_mismatch]
    
    if (length(mismatch_positions) == 0) {
      return(list())
    }
    
    # Map alignment positions back to donor sequence positions
    donor_aligned <- strsplit(attr(mismatch_detail, "alignment")$aligned_ref, "")[[1]]
    
    # Create mapping from alignment position to sequence position
    seq_pos <- 0
    align_to_seq <- integer(length(donor_aligned))
    for (i in seq_along(donor_aligned)) {
      if (donor_aligned[i] != "-") {
        seq_pos <- seq_pos + 1
      }
      align_to_seq[i] <- seq_pos
    }
    
    # Convert alignment positions to sequence positions
    mismatch_seq_positions <- align_to_seq[mismatch_positions]
    mismatch_seq_positions <- mismatch_seq_positions[mismatch_seq_positions > 0]
    
    # Remove gaps from donor sequence for peptide generation
    donor_seq_nogaps <- gsub("-", "", attr(mismatch_detail, "alignment")$aligned_ref)
    seq_len <- nchar(donor_seq_nogaps)
    
    peptide_data <- list()
    
    # Generate peptides around each mismatch position for ALL lengths
    all_lengths <- c(class_I_lengths, class_II_lengths)
    
    for (pos in unique(mismatch_seq_positions)) {
      for (len in all_lengths) {
        # Determine MHC class for this length
        mhc_class <- if (len %in% class_I_lengths) "I" else "II"
        
        # Calculate window of starts that would include this position
        min_start <- max(1, pos - len + 1)
        max_start <- min(seq_len - len + 1, pos)
        
        if (max_start >= min_start) {
          for (start in min_start:max_start) {
            pep <- substr(donor_seq_nogaps, start, start + len - 1)
            if (nchar(pep) == len) {
              peptide_data[[length(peptide_data) + 1]] <- list(
                peptide = pep,
                start = start,
                end = start + len - 1,
                length = len,
                mhc_class = mhc_class,
                mismatch_position = pos,
                position_in_peptide = pos - start + 1,
                compared_to = recipient_allele,
                donor_allele = donor_allele
              )
            }
          }
        }
      }
    }
    
    peptide_data
  }
  
  # Helper function to find most similar recipient allele at same locus
  findClosestRecipient <- function(donor_allele, donor_seq, recipient_seqs, recipient_alleles) {
    # Get locus of donor allele
    donor_locus <- sub("\\*.*", "", donor_allele)
    
    # Find recipient alleles at same locus
    same_locus <- grep(paste0("^", donor_locus, "\\*"), recipient_alleles)
    
    if (length(same_locus) == 0) {
      # No recipient alleles at same locus - use all
      return(recipient_alleles)
    }
    
    if (length(same_locus) == 1) {
      return(recipient_alleles[same_locus[1]])
    }
    
    # Compare to each recipient at same locus and find the one with fewest mismatches
    mismatch_counts <- sapply(same_locus, function(i) {
      quantifyMismatch(donor_seq, recipient_seqs[[i]], return = "count")
    })
    
    recipient_alleles[same_locus[which.min(mismatch_counts)]]
  }
  
  # Step 4: Generate ALL peptides from ALL mismatched alleles
  message("=== Generating peptides from mismatched alleles ===")
  
  all_peptide_data <- list()
  comparison_summary <- list()
  
  for (donor_allele in mismatched_alleles) {
    donor_seq <- prot_seq_map[[donor_allele]]
    
    if (is.null(donor_seq) || nchar(donor_seq) == 0) {
      warning(sprintf("Could not retrieve sequence for %s", donor_allele))
      next
    }
    
    # Get locus of donor allele
    donor_locus <- sub("\\*.*", "", donor_allele)
    
    # Determine which recipient allele(s) to compare against for peptide generation
    # (This is just for alignment to find mismatches - we'll test against ALL later)
    if (compare_to == "closest") {
      # Find recipients at same locus
      locus_pattern <- paste0("^", donor_locus, "\\*")
      recipient_locus <- r_alleles[grep(locus_pattern, r_alleles)]
      
      if (length(recipient_locus) > 0) {
        recipient_seqs <- lapply(recipient_locus, function(a) prot_seq_map[[a]])
        names(recipient_seqs) <- recipient_locus
        compare_alleles <- findClosestRecipient(donor_allele, donor_seq, recipient_seqs, recipient_locus)
      } else {
        # No recipients at same locus - use first available
        compare_alleles <- r_alleles[1]
      }
    } else {
      # Compare to all recipients at same locus
      locus_pattern <- paste0("^", donor_locus, "\\*")
      compare_alleles <- r_alleles[grep(locus_pattern, r_alleles)]
      
      if (length(compare_alleles) == 0) {
        # No recipients at same locus - use all
        compare_alleles <- r_alleles
      }
    }
    
    # Generate peptides by comparing to selected recipient allele(s)
    for (compare_allele in compare_alleles) {
      recipient_seq <- prot_seq_map[[compare_allele]]
      
      peptide_data <- generatePeptidesFromMismatches(
        donor_seq, recipient_seq,
        donor_allele, compare_allele,
        class_I_lengths = 8:11,
        class_II_lengths = 12:25
      )
      
      if (length(peptide_data) > 0) {
        all_peptide_data <- c(all_peptide_data, peptide_data)
        
        comparison_summary[[paste(donor_allele, compare_allele, sep = "_vs_")]] <- list(
          donor = donor_allele,
          recipient = compare_allele,
          n_mismatches = length(unique(sapply(peptide_data, function(x) x$mismatch_position))),
          n_peptides = length(peptide_data)
        )
      }
    }
    
    message(sprintf("  %s: Generated peptides from alignment comparisons", donor_allele))
  }
  
  if (length(all_peptide_data) == 0) {
    stop("No peptides generated from mismatched alleles.")
  }
  
  # Create peptide metadata dataframe
  message("=== Creating peptide database ===")
  peptide_metadata <- do.call(rbind, lapply(all_peptide_data, function(pd) {
    data.frame(
      peptide = pd$peptide,
      peptide_length = pd$length,
      mhc_class = pd$mhc_class,
      start_position = pd$start,
      end_position = pd$end,
      mismatch_position = pd$mismatch_position,
      position_in_peptide = pd$position_in_peptide,
      compared_to = pd$compared_to,
      donor_allele = pd$donor_allele,
      donor_locus = sub("\\*.*", "", pd$donor_allele),
      stringsAsFactors = FALSE
    )
  }))
  
  # Get unique peptides by class
  class_I_peptides <- unique(peptide_metadata$peptide[peptide_metadata$mhc_class == "I"])
  class_II_peptides <- unique(peptide_metadata$peptide[peptide_metadata$mhc_class == "II"])
  
  message(sprintf("Generated %d unique Class I peptides (8-11mers)", length(class_I_peptides)))
  message(sprintf("Generated %d unique Class II peptides (12-25mers)", length(class_II_peptides)))
  
  # Step 5: Test ALL peptides against ALL recipient alleles (with appropriate class matching)
  message("=== Testing peptides against recipient alleles ===")
  
  all_results <- list()
  
  # Test Class I peptides against Class I recipient alleles
  if (length(class_I_peptides) > 0 && length(class_I_recipients) > 0) {
    message(sprintf("\nTesting %d Class I peptides against %d Class I recipient alleles",
                    length(class_I_peptides), length(class_I_recipients)))
    
    for (recipient_allele in class_I_recipients) {
      tryCatch({
        predictions <- predictMHCnuggets(
          peptides = class_I_peptides,
          allele = recipient_allele,
          mhc_class = "I",
          ic50_threshold = ic50_threshold,
          rank_output = TRUE,
          hla_env = hla_env
        )
        
        # Merge with metadata (keeping all metadata)
        predictions <- merge(predictions, 
                             peptide_metadata[peptide_metadata$mhc_class == "I", ], 
                             by = "peptide", all.x = TRUE)
        
        predictions$recipient_allele <- recipient_allele
        predictions$recipient_locus <- sub("\\*.*", "", recipient_allele)
        predictions$binding <- predictions$ic50 <= ic50_threshold
        
        all_results[[paste("ClassI", recipient_allele, sep = "_")]] <- predictions
        
        message(sprintf("  %s: %d/%d peptides bind (IC50 <= %d nM)",
                        recipient_allele, 
                        sum(predictions$binding),
                        nrow(predictions),
                        ic50_threshold))
        
      }, error = function(e) {
        warning(sprintf("Failed to predict for allele %s: %s", recipient_allele, e$message))
      })
    }
  }
  
  # Test Class II peptides against Class II recipient alleles
  if (length(class_II_peptides) > 0 && length(class_II_recipients) > 0) {
    message(sprintf("\nTesting %d Class II peptides against %d Class II recipient alleles",
                    length(class_II_peptides), length(class_II_recipients)))
    
    class_II_pairs <- get_classII_pairs(class_II_recipients)
    for (recipient_allele in class_II_pairs) {
      tryCatch({
        predictions <- predictMHCnuggets(
          peptides = class_II_peptides,
          allele = recipient_allele,
          mhc_class = "II",
          ic50_threshold = ic50_threshold,
          rank_output = TRUE,
          hla_env = hla_env
        )
        
        # Merge with metadata (keeping all metadata)
        predictions <- merge(predictions, 
                             peptide_metadata[peptide_metadata$mhc_class == "II", ], 
                             by = "peptide", all.x = TRUE)
        
        predictions$recipient_allele <- recipient_allele
        predictions$recipient_locus <- sub("\\*.*", "", recipient_allele)
        predictions$binding <- predictions$ic50 <= ic50_threshold
        
        all_results[[paste("ClassII", recipient_allele, sep = "_")]] <- predictions
        
        message(sprintf("  %s: %d/%d peptides bind (IC50 <= %d nM)",
                        recipient_allele, 
                        sum(predictions$binding),
                        nrow(predictions),
                        ic50_threshold))
        
      }, error = function(e) {
        warning(sprintf("Failed to predict for allele %s: %s", recipient_allele, e$message))
      })
    }
  }
  
  if (length(all_results) == 0) {
    stop("No predictions were generated.")
  }
  
  # Combine all results
  all_predictions <- do.call(rbind, all_results)
  
  # Format return value based on request
  if (return == "summary") {
    total_peptides <- nrow(all_predictions)
    binding_peptides <- sum(all_predictions$binding)
    
    summary_df <- data.frame(
      total_predictions = total_peptides,
      unique_peptides_tested = length(unique(all_predictions$peptide)),
      binding_predictions = binding_peptides,
      binding_percentage = if (total_peptides > 0) 100 * binding_peptides / total_peptides else 0,
      ic50_threshold = ic50_threshold,
      n_mismatched_alleles = length(mismatched_alleles),
      n_recipient_alleles = length(r_alleles),
      n_class_I_recipients = length(class_I_recipients),
      n_class_II_recipients = length(class_II_recipients)
    )
    
    return(summary_df)
    
  } else if (return == "by_recipient_allele") {
    # Group by recipient allele showing total bound peptides
    by_recipient <- all_predictions |>
      group_by(recipient_allele, recipient_locus, mhc_class) |>
      summarise(
        n_peptides_tested = n(),
        n_unique_peptides = n_distinct(peptide),
        n_binding_peptides = sum(binding),
        binding_percentage = 100 * n_binding_peptides / n_peptides_tested,
        mean_ic50_all = mean(ic50, na.rm = TRUE),
        mean_ic50_binders = mean(ic50[binding], na.rm = TRUE),
        median_ic50_binders = median(ic50[binding], na.rm = TRUE),
        n_donor_alleles = n_distinct(donor_allele),
        donor_loci = paste(unique(donor_locus), collapse = ", "),
        .groups = "drop"
      ) |>
      arrange(desc(n_binding_peptides))
    
    return(as.data.frame(by_recipient))
    
  } else {  # return == "detailed"
    return(list(
      summary = data.frame(
        total_predictions = nrow(all_predictions),
        unique_peptides_tested = length(unique(all_predictions$peptide)),
        binding_predictions = sum(all_predictions$binding),
        n_mismatched_alleles = length(mismatched_alleles),
        n_recipient_alleles = length(r_alleles)
      ),
      mismatched_alleles = mismatched_alleles,
      recipient_alleles = r_alleles,
      class_I_recipients = class_I_recipients,
      class_II_recipients = class_II_recipients,
      comparison_summary = comparison_summary,
      all_predictions = all_predictions
    ))
  }
}


#' Visualize Cross-Locus Peptide Binding Results
#'
#' @description
#' Creates visualizations of peptide binding predictions across all loci
#'
#' @param binding_results Results from calculatePeptideBindingLoad with return="detailed"
#' @param plot_type Type of plot: "heatmap", "bar_by_recipient", "bar_by_donor", or "scatter"
#' @param palette Character. A color palette name. Defaults to "spectral".
#' @param ... Additional arguments passed to the ggplot theme.
#'
#' @return ggplot object
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_bar geom_point scale_fill_gradient2 theme_minimal labs
#' @importFrom dplyr group_by summarise
#' @export
visualizePeptideBinding <- function(binding_results, 
                                    plot_type = c("heatmap", "bar_by_recipient", "bar_by_donor", "scatter"), 
                                    palette = "spectral", 
                                    ...) {
  plot_type <- match.arg(plot_type)
  
  if (!is.list(binding_results) || !"all_predictions" %in% names(binding_results)) {
    stop("binding_results must be output from calculatePeptideBindingLoad with return='detailed'")
  }
  
  data <- binding_results$all_predictions
  
  if (plot_type == "heatmap") {
    # Cross-locus heatmap: donor alleles vs recipient alleles
    summary_data <- data |>
      dplyr::group_by(donor_allele, recipient_allele) |>
      dplyr::summarise(
        binding_rate = mean(binding) * 100,
        n_peptides = dplyr::n(),
        .groups = "drop"
      )
    
    p <- ggplot2::ggplot(summary_data, ggplot2::aes(x = donor_allele, y = recipient_allele, fill = binding_rate)) +
      ggplot2::geom_tile(color = "white", lwd = 0.5) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%\n(%d)", binding_rate, n_peptides)),
                         size = 3, color = "black") +
      ggplot2::scale_fill_gradientn(colors = rev(.colorizer(n=11, palette = palette)),
                                    name = "Binding %",
                                    limits = c(0, 100)) +
      .themeMatchR(...) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 9)
      ) +
      ggplot2::labs(title = "Cross-Locus Peptide Binding: All Mismatched Peptides vs All Recipients",
                    subtitle = "Percentage of binding peptides (number tested)",
                    x = "Mismatched Donor Allele", 
                    y = "Recipient Allele")
    
  } else if (plot_type == "bar_by_recipient") {
    # Bar plot by recipient allele showing total bound peptides
    recipient_summary <- data |>
      dplyr::group_by(recipient_allele, recipient_locus, mhc_class) |>
      dplyr::summarise(
        total = dplyr::n(),
        binding = sum(binding),
        .groups = "drop"
      ) |>
      dplyr::arrange(desc(binding))
    
    p <- ggplot2::ggplot(recipient_summary, ggplot2::aes(x = reorder(recipient_allele, binding), 
                                                         y = binding, 
                                                         fill = mhc_class)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::geom_text(ggplot2::aes(label = paste0(binding, "/", total)), 
                         hjust = -0.1, size = 3) +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("I" = .colorizer(n=2, palette = palette)[1],
                                            "II" = .colorizer(n=2, palette = palette)[2]),
                                 name = "MHC Class") +
      .themeMatchR(...) +
      ggplot2::labs(title = "Binding Peptides by Recipient Allele",
                    subtitle = "Total bound mismatched peptides from all donor alleles",
                    x = "Recipient Allele", 
                    y = "Number of Binding Peptides")
    
  } else if (plot_type == "bar_by_donor") {
    # Bar plot by donor allele
    donor_summary <- data |>
      dplyr::group_by(donor_allele, donor_locus) |>
      dplyr::summarise(
        total = dplyr::n(),
        binding = sum(binding),
        n_recipient_alleles = dplyr::n_distinct(recipient_allele),
        .groups = "drop"
      ) |>
      dplyr::arrange(desc(binding))
    
    p <- ggplot2::ggplot(donor_summary, ggplot2::aes(x = reorder(donor_allele, binding), 
                                                     y = binding)) +
      ggplot2::geom_bar(stat = "identity", fill = .colorizer(n=2, palette = palette)[2]) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%d/%d\n(%d alleles)", 
                                                      binding, total, n_recipient_alleles)), 
                         hjust = -0.1, size = 3) +
      ggplot2::coord_flip() +
      .themeMatchR(...) +
      ggplot2::labs(title = "Binding Peptides by Donor Allele",
                    subtitle = "Total peptides binding to any recipient allele",
                    x = "Mismatched Donor Allele", 
                    y = "Number of Binding Peptides")
    
  } else {  # scatter
    # Scatter plot of IC50 values by recipient allele
    p <- ggplot2::ggplot(data, ggplot2::aes(x = ic50, y = recipient_allele, color = binding)) +
      ggplot2::geom_point(alpha = 0.4, position = ggplot2::position_jitter(height = 0.2)) +
      ggplot2::scale_x_log10() +
      ggplot2::geom_vline(xintercept = 500, linetype = "dashed", 
                          color = .colorizer(n=2, palette = palette)[1]) +
      ggplot2::scale_color_manual(values = c("FALSE" = "gray", 
                                             "TRUE" = .colorizer(n=2, palette = palette)[1]),
                                  name = "Binding") +
      ggplot2::facet_wrap(~mhc_class, scales = "free_y") +
      .themeMatchR(...) +
      ggplot2::labs(title = "IC50 Distribution by Recipient Allele",
                    subtitle = "All mismatched peptides tested",
                    x = "IC50 (nM, log scale)", 
                    y = "Recipient Allele")
  }
  
  return(p)
}

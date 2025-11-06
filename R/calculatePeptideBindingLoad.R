#' Calculate Peptide Binding Load from Mismatched Positions
#'
#' @description
#' Generates peptides from mismatched positions between donor and recipient
#' and predicts MHC binding using MHCnuggets. Comprehensively examines all
#' mismatched peptides between every mismatched allele pair.
#' 
#' @param recipient_geno Recipient HLA genotype object
#' @param donor_geno Donor HLA genotype object
#' @param loci Character vector of loci to analyze (NULL for all)
#' @param peptide_lengths For Class I: integer vector of peptide lengths (default 8:11)
#' @param mhc_class "I" or "II" - determines peptide length defaults
#' @param ic50_threshold IC50 threshold in nM for binding (default 500)
#' @param return What to return: "summary", "detailed", or "per_locus"
#' @param parallel Use parallel processing
#' @param n_cores Number of cores for parallel processing
#' @param hla_env Basilisk environment for MHCnuggets
#' @param cache_dir Optional directory for sequence caching
#'
#' @return Enhanced results including nucleotide analysis if requested
#'
#' @export
calculatePeptideBindingLoad <- function(recipient_geno,
                                        donor_geno,
                                        loci = NULL,
                                        peptide_lengths = NULL,
                                        mhc_class = "I",
                                        ic50_threshold = 500,
                                        return = c("summary", "detailed", "per_locus"),
                                        parallel = TRUE,
                                        n_cores = NULL,
                                        hla_env = deepmatchrEnv(),
                                        cache_dir = NULL) {
  
  return <- match.arg(return)
  mhc_class <- toupper(mhc_class)
  
  # Set default peptide lengths based on MHC class
  if (is.null(peptide_lengths)) {
    peptide_lengths <- if (mhc_class == "I") 8:11 else 12:25
  }
  
  # Validate genotypes
  validateHlaGeno(recipient_geno)
  validateHlaGeno(donor_geno)
  
  # Determine loci to analyze
  shared_loci <- intersect(recipient_geno$locus_present, donor_geno$locus_present)
  if (!is.null(loci)) shared_loci <- intersect(shared_loci, loci)
  if (length(shared_loci) == 0) stop("No shared loci between donor and recipient.")
  
  # Extract alleles
  get_alleles_for_loci <- function(geno, loci_vec) {
    pat <- paste0("^(", paste(loci_vec, collapse = "|"), ")_")
    idx <- grep(pat, colnames(geno$data))
    alleles <- unlist(geno$data[1, idx, drop = FALSE])
    alleles[!is.na(alleles) & nzchar(alleles)]
  }
  
  r_alleles <- get_alleles_for_loci(recipient_geno, shared_loci)
  d_alleles <- get_alleles_for_loci(donor_geno, shared_loci)
  
  # Batch retrieve protein sequences
  all_alleles <- unique(c(r_alleles, d_alleles))
  prot_seq_map <- batchGetSequences(all_alleles, type = "PROT", 
                                    n_cores = n_cores, cache_dir = cache_dir)
  
  # Enhanced mismatch analysis function
  analyzeMismatches <- function(seq1_prot, seq2_prot) {
    s1 <- unlist(strsplit(seq1_prot, ""))
    s2 <- unlist(strsplit(seq2_prot, ""))
    max_len <- max(length(s1), length(s2))
    
    # Pad shorter sequence
    if (length(s1) < max_len) s1 <- c(s1, rep("-", max_len - length(s1)))
    if (length(s2) < max_len) s2 <- c(s2, rep("-", max_len - length(s2)))
    
    # Find positions where sequences differ
    mismatch_positions <- which(s1 != s2)
    
    # If nucleotide sequences provided, analyze codon changes
    mismatch_details <- data.frame(
      position = mismatch_positions,
      aa_donor = s1[mismatch_positions],
      aa_recipient = s2[mismatch_positions],
      stringsAsFactors = FALSE
    )
    
    mismatch_details
  }
  
  # Function to generate peptides with position tracking
  generatePeptidesWithPositions <- function(sequence, position, lengths) {
    seq_len <- nchar(sequence)
    peptide_data <- list()
    
    for (len in lengths) {
      min_start <- max(1, position - len + 1)
      max_start <- min(seq_len - len + 1, position)
      
      if (max_start >= min_start) {
        for (start in min_start:max_start) {
          pep <- substr(sequence, start, start + len - 1)
          if (nchar(pep) == len) {
            peptide_data[[length(peptide_data) + 1]] <- list(
              peptide = pep,
              start = start,
              end = start + len - 1,
              mismatch_position = position,
              contains_position = position - start + 1  # Position within peptide
            )
          }
        }
      }
    }
    
    peptide_data
  }
  
  # Process each locus with enhanced analysis
  process_locus_enhanced <- function(L) {
    rL <- r_alleles[grep(paste0("^", L, "_"), names(r_alleles))]
    dL <- d_alleles[grep(paste0("^", L, "_"), names(d_alleles))]
    
    if (length(rL) == 0 || length(dL) == 0) {
      return(list(
        locus = L,
        total_peptides = 0,
        binding_peptides = 0,
        results = data.frame(),
        mismatch_analysis = data.frame()
      ))
    }
    
    all_results <- list()
    all_mismatches <- list()
    
    # Compare each recipient-donor allele pair
    for (r_allele in rL) {
      for (d_allele in dL) {
        if (r_allele == d_allele) next
        
        r_seq_prot <- prot_seq_map[[r_allele]]
        d_seq_prot <- prot_seq_map[[d_allele]]
        
        
        # Analyze mismatches
        mismatch_analysis <- analyzeMismatches(
          d_seq_prot, r_seq_prot, d_seq_nuc, r_seq_nuc
        )
        
        if (nrow(mismatch_analysis) == 0) next
        
        mismatch_analysis$donor_allele <- d_allele
        mismatch_analysis$recipient_allele <- r_allele
        mismatch_analysis$locus <- L
        
        all_mismatches[[paste(r_allele, d_allele, sep = "_")]] <- mismatch_analysis
        
        # Generate peptides with position tracking
        donor_peptide_data <- unlist(lapply(mismatch_analysis$position, function(pos) {
          generatePeptidesWithPositions(d_seq_prot, pos, peptide_lengths)
        }), recursive = FALSE)
        
        recipient_peptide_data <- unlist(lapply(mismatch_analysis$position, function(pos) {
          generatePeptidesWithPositions(r_seq_prot, pos, peptide_lengths)
        }), recursive = FALSE)
        
        # Extract unique donor peptides
        donor_peptides <- unique(sapply(donor_peptide_data, function(x) x$peptide))
        recipient_peptides <- unique(sapply(recipient_peptide_data, function(x) x$peptide))
        unique_donor_peptides <- setdiff(donor_peptides, recipient_peptides)
        
        if (length(unique_donor_peptides) > 0) {
          # Create peptide metadata
          peptide_metadata <- do.call(rbind, lapply(donor_peptide_data, function(pd) {
            if (pd$peptide %in% unique_donor_peptides) {
              data.frame(
                peptide = pd$peptide,
                start_position = pd$start,
                end_position = pd$end,
                mismatch_position = pd$mismatch_position,
                mismatch_position_in_peptide = pd$contains_position,
                stringsAsFactors = FALSE
              )
            }
          }))
          
          # Remove NULL entries and duplicates
          if (!is.null(peptide_metadata)) {
            peptide_metadata <- unique(peptide_metadata)
            
            # Predict binding for recipient alleles
            for (test_allele in rL) {
              tryCatch({
                predictions <- predictMHCnuggets(
                  peptides = unique_donor_peptides,
                  allele = test_allele,
                  mhc_class = mhc_class,
                  ic50_threshold = ic50_threshold,
                  rank_output = TRUE,
                  hla_env = hla_env
                )
                
                # Merge with metadata
                predictions <- merge(predictions, peptide_metadata, 
                                     by = "peptide", all.x = TRUE)
                
                predictions$donor_allele <- d_allele
                predictions$recipient_allele <- r_allele
                predictions$test_allele <- test_allele
                predictions$locus <- L
                predictions$binding <- predictions$ic50 <= ic50_threshold
                
                
                all_results[[paste(r_allele, d_allele, test_allele, sep = "_")]] <- predictions
                
              }, error = function(e) {
                warning(sprintf("Failed to predict for allele %s: %s", test_allele, e$message))
              })
            }
          }
        }
      }
    }
    
    # Combine results
    if (length(all_results) > 0) {
      locus_results <- do.call(rbind, all_results)
      mismatch_summary <- if (length(all_mismatches) > 0) {
        do.call(rbind, all_mismatches)
      } else {
        data.frame()
      }
      
      return(list(
        locus = L,
        total_peptides = nrow(locus_results),
        binding_peptides = sum(locus_results$binding),
        mean_ic50_binders = mean(locus_results$ic50[locus_results$binding]),
        results = locus_results,
        mismatch_analysis = mismatch_summary
      ))
    } else {
      return(list(
        locus = L,
        total_peptides = 0,
        binding_peptides = 0,
        mean_ic50_binders = NA,
        results = data.frame(),
        mismatch_analysis = if (length(all_mismatches) > 0) do.call(rbind, all_mismatches) else data.frame()
      ))
    }
  }
  
  # Process all loci
  if (parallel && length(shared_loci) > 1 && !is.null(n_cores) && n_cores > 1) {
    locus_results <- parallel::mclapply(shared_loci, process_locus_enhanced, mc.cores = n_cores)
  } else {
    locus_results <- lapply(shared_loci, process_locus_enhanced)
  }
  
  # Format return value based on request
  if (return == "summary") {
    total_peptides <- sum(sapply(locus_results, function(x) x$total_peptides))
    binding_peptides <- sum(sapply(locus_results, function(x) x$binding_peptides))
    
    summary_df <- data.frame(
      total_mismatched_peptides = total_peptides,
      binding_peptides = binding_peptides,
      binding_percentage = if (total_peptides > 0) 100 * binding_peptides / total_peptides else 0,
      ic50_threshold = ic50_threshold,
      mhc_class = mhc_class,
    )
    
    return(summary_df)
    
  } else if (return == "per_locus") {
    return(do.call(rbind, lapply(locus_results, function(x) {
      df <- data.frame(
        locus = x$locus,
        total_peptides = x$total_peptides,
        binding_peptides = x$binding_peptides,
        binding_percentage = if (x$total_peptides > 0) 100 * x$binding_peptides / x$total_peptides else 0,
        mean_ic50_binders = x$mean_ic50_binders
      )
      
      df
    })))
    
  } else {  # return == "detailed"
    return(list(
      summary = data.frame(
        total_mismatched_peptides = sum(sapply(locus_results, function(x) x$total_peptides)),
        binding_peptides = sum(sapply(locus_results, function(x) x$binding_peptides))
      ),
      per_locus = locus_results,
      all_predictions = do.call(rbind, lapply(locus_results, function(x) x$results)),
    ))
  }
}

#' Visualize Peptide Binding Results
#'
#' @description
#' Creates visualizations of peptide binding predictions
#'
#' @param binding_results Results from calculatePeptideBindingLoad with return="detailed"
#' @param plot_type Type of plot: "heatmap", "bar", or "scatter"
#'
#' @return ggplot object
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_bar geom_point scale_fill_gradient2 theme_minimal labs
#' @importFrom dplyr group_by summarise
#' @export
visualizePeptideBinding <- function(binding_results, plot_type = c("heatmap", "bar", "scatter")) {
  plot_type <- match.arg(plot_type)
  
  if (!is.list(binding_results) || !"all_predictions" %in% names(binding_results)) {
    stop("binding_results must be output from calculatePeptideBindingLoad with return='detailed'")
  }
  
  data <- binding_results$all_predictions
  
  if (plot_type == "heatmap") {
    # Aggregate by allele pairs
    summary_data <- data %>%
      dplyr::group_by(donor_allele, test_allele, locus) %>%
      dplyr::summarise(
        binding_rate = mean(binding) * 100,
        mean_ic50 = mean(ic50[binding]),
        .groups = "drop"
      )
    
    p <- ggplot2::ggplot(summary_data, ggplot2::aes(x = donor_allele, y = test_allele, fill = binding_rate)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 50,
                                    name = "Binding %") +
      ggplot2::facet_wrap(~locus, scales = "free") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = "Peptide Binding Rates by Allele Pair",
                    x = "Donor Allele", y = "Recipient Test Allele")
    
  } else if (plot_type == "bar") {
    # Bar plot by locus
    locus_summary <- data %>%
      dplyr::group_by(locus) %>%
      dplyr::summarise(
        total = n(),
        binding = sum(binding),
        .groups = "drop"
      )
    
    p <- ggplot2::ggplot(locus_summary, ggplot2::aes(x = locus, y = binding)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::geom_text(ggplot2::aes(label = paste0(binding, "/", total)), 
                         vjust = -0.5) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Binding Peptides by Locus",
                    x = "Locus", y = "Number of Binding Peptides")
    
  } else {  # scatter
    # Scatter plot of IC50 values
    p <- ggplot2::ggplot(data, ggplot2::aes(x = ic50, y = test_allele, color = binding)) +
      ggplot2::geom_point(alpha = 0.6, position = ggplot2::position_jitter(height = 0.2)) +
      ggplot2::scale_x_log10() +
      ggplot2::geom_vline(xintercept = 500, linetype = "dashed", color = "red") +
      ggplot2::scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                                  name = "Binding") +
      ggplot2::facet_wrap(~locus, scales = "free_y") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "IC50 Distribution by Allele",
                    x = "IC50 (nM, log scale)", y = "Test Allele")
  }
  
  return(p)
}
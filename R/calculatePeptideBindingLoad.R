#' Calculate Peptide Binding Load for Transplant Risk Assessment
#'
#' @description
#' Predicts transplant risk by calculating peptide-HLA binding affinities
#' between recipient HLA molecules and donor-mismatched peptides. Supports
#' multiple binding prediction backends: built-in position weight matrix (PWM),
#' NetMHCpan, or MHCnuggets.
#'
#' @param recipient An `hla_genotype` object or character vector of HLA allele names.
#' @param donor An `hla_genotype` object, character vector of HLA allele names, or
#'   a character vector of peptide sequences. If `hla_genotype` or allele names,
#'   mismatched peptides are derived automatically from sequence differences.
#' @param backend Character. Binding prediction method: `"pwm"` (default, no external
#'   dependencies), `"netmhcpan"`, or `"mhcnuggets"`.
#' @param backend_path Character. Path to external tool executable. Required for
#'
#'   `"netmhcpan"` backend. Download NetMHCpan from
#'   \url{https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/}.
#' @param peptide_length Integer. Peptide length(s) to consider. Default `9L`.
#' @param binding_threshold Numeric. IC50 threshold (nM) for "strong binder".
#'   Default `500`.
#' @param weak_threshold Numeric. IC50 threshold (nM) for "weak binder".
#'   Default `5000`.
#' @param return Character. One of `"total"` (risk score), `"summary"` (per-allele),
#'   or `"detail"` (per-peptide table). Default `"total"`.
#' @param aggregate_method Character. How to combine per-peptide scores:
#'   `"sum"`, `"max"`, `"mean"`. Default `"sum"`.
#'
#' @return Depends on `return`:
#'   - `"total"`: Numeric risk score.
#'   - `"summary"`: data.frame with per-HLA-allele binding summary.
#'   - `"detail"`: data.frame with columns: `peptide`, `hla_allele`, `predicted_ic50`,
#'     `binding_level`, `contribution`.
#'
#' @details
#' The function works in several steps:
#' 1. **Input processing**: Converts inputs to standard format (allele names and peptides)
#' 2. **Peptide derivation**: If donor is genotype/alleles, derives mismatched peptides
#'    by comparing sequences and generating overlapping k-mers from mismatch regions

#' 3. **Binding prediction**: Uses selected backend to predict IC50 values
#' 4. **Risk calculation**: Aggregates binding predictions into a risk score
#'
#' The **PWM backend** uses simplified position weight matrices based on HLA supertypes.
#' For production use with high accuracy requirements, NetMHCpan or MHCnuggets is recommended.
#'
#' **External backends:**
#' - **NetMHCpan**: A state-of-the-art method for predicting peptide-MHC class I binding
#'   using artificial neural networks. Available at
#'   \url{https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/}.
#' - **MHCnuggets**: A deep learning approach for MHC binding prediction. Available at
#'   \url{https://github.com/KarchinLab/mhcnuggets}. See \code{\link{predictMHCnuggets}}
#'   for direct access to MHCnuggets predictions.
#'
#' Risk score formula:
#' \deqn{contribution = (1 - IC50/weak\_threshold) \times multiplier}
#' where multiplier is 2 for strong binders, 1 for weak binders.
#'
#' @examples
#' # Example 1: Using raw peptides (no external data required)
#' # Define recipient HLA alleles
#' recipient_alleles <- c("A*02:01", "A*03:01", "B*07:02", "B*08:01")
#'
#' # Define peptides to test
#' peptides <- c("GILGFVFTL", "NLVPMVATV", "FLKEKGGL", "SIINFEKL")
#'
#' # Calculate binding load with PWM backend
#' result <- calculatePeptideBindingLoad(
#'   recipient = recipient_alleles,
#'   donor = peptides,
#'   backend = "pwm",
#'   return = "summary"
#' )
#' print(result)
#'
#' # Get detailed per-peptide results
#' detail <- calculatePeptideBindingLoad(
#'   recipient = recipient_alleles,
#'   donor = peptides,
#'   backend = "pwm",
#'   return = "detail"
#' )
#' head(detail)
#'
#' # Example 2: Using hla_genotype objects with peptides
#' recipient <- data.frame(
#'   A_1 = "A*02:01", A_2 = "A*03:01",
#'   B_1 = "B*07:02", B_2 = "B*44:02"
#' )
#' rgeno <- hlaGeno(recipient)
#'
#' # Calculate total risk score
#' total_risk <- calculatePeptideBindingLoad(
#'   recipient = rgeno,
#'   donor = peptides,
#'   return = "total"
#' )
#' print(total_risk)
#'
#' \donttest{
#' # Example 3: Using genotypes (requires IMGT database connection)
#' donor <- data.frame(
#'   A_1 = "A*01:01", A_2 = "A*24:02",
#'   B_1 = "B*08:01", B_2 = "B*35:01"
#' )
#' dgeno <- hlaGeno(donor)
#'
#' # Calculate binding load from sequence mismatches
#' calculatePeptideBindingLoad(rgeno, dgeno, return = "summary")
#' }
#'
#' @references
#' Reynisson B, et al. (2020). NetMHCpan-4.1 and NetMHCIIpan-4.0: improved 
#' predictions of MHC antigen presentation by concurrent motif deconvolution 
#' and integration of MS MHC eluted ligand data. *Nucleic Acids Research*, 
#' 48(W1), W449-W454. \doi{10.1093/nar/gkaa379}
#'
#' Shao XM, et al. (2020). High-Throughput Prediction of MHC Class I and II
#' Neoantigens with MHCnuggets. *Cancer Immunology Research*, 8(3), 396-408.
#' \doi{10.1158/2326-6066.CIR-19-0464}
#'
#' @seealso \code{\link{calculateMismatchLoad}}, \code{\link{quantifyMismatch}},
#'   \code{\link{predictMHCnuggets}}
#'
#' @export
calculatePeptideBindingLoad <- function(
    recipient,
    donor,
    backend = c("pwm", "netmhcpan", "mhcnuggets"),
    backend_path = NULL,
    peptide_length = 9L,
    binding_threshold = 500,
    weak_threshold = 5000,
    return = c("total", "summary", "detail"),
    aggregate_method = c("sum", "max", "mean")
) {
  backend <- match.arg(backend)
  return <- match.arg(return)
  aggregate_method <- match.arg(aggregate_method)

 # --- 1. Process recipient to get HLA alleles ---
  recipient_alleles <- .extractAlleles(recipient)
  if (length(recipient_alleles) == 0L) {
    stop("No valid recipient HLA alleles found.")
  }

  # --- 2. Process donor to get peptides ---
  peptides <- .getPeptides(donor, recipient, peptide_length)
  if (length(peptides) == 0L) {
    if (return == "total") return(0)
    if (return == "summary") {
      return(data.frame(
        hla_allele = recipient_alleles,
        n_peptides = 0L,
        n_strong = 0L,
        n_weak = 0L,
        risk_contribution = 0,
        stringsAsFactors = FALSE
      ))
    }
    return(data.frame(
      peptide = character(0),
      hla_allele = character(0),
      predicted_ic50 = numeric(0),
      binding_level = character(0),
      contribution = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  # --- 3. Predict binding for each peptide x allele combination ---
  results <- .predictBinding(
    peptides = peptides,
    alleles = recipient_alleles,
    backend = backend,
    backend_path = backend_path,
    peptide_length = peptide_length
  )

  # --- 4. Calculate binding levels and contributions ---
  results$binding_level <- ifelse(
    results$predicted_ic50 <= binding_threshold, "strong",
    ifelse(results$predicted_ic50 <= weak_threshold, "weak", "non_binder")
  )

  results$contribution <- ifelse(
    results$binding_level == "strong",
    2 * (1 - results$predicted_ic50 / weak_threshold),
    ifelse(
      results$binding_level == "weak",
      1 * (1 - results$predicted_ic50 / weak_threshold),
      0
    )
  )
  results$contribution <- pmax(results$contribution, 0)

  # --- 5. Return based on requested format ---
  if (return == "detail") {
    return(results)
  }

  if (return == "summary") {
    summary_df <- do.call(rbind, lapply(recipient_alleles, function(a) {
      subset <- results[results$hla_allele == a, ]
      data.frame(
        hla_allele = a,
        n_peptides = nrow(subset),
        n_strong = sum(subset$binding_level == "strong"),
        n_weak = sum(subset$binding_level == "weak"),
        risk_contribution = sum(subset$contribution),
        stringsAsFactors = FALSE
      )
    }))
    return(summary_df)
  }

  # return == "total"
  total_risk <- switch(aggregate_method,
    sum = sum(results$contribution),
    max = if (nrow(results) > 0) max(results$contribution) else 0,
    mean = if (nrow(results) > 0) mean(results$contribution) else 0
  )
  return(total_risk)
}


# --- Internal helper functions ---

#' Extract allele names from various input types
#' @keywords internal
.extractAlleles <- function(x) {
  if (inherits(x, "hla_genotype")) {
    validateHlaGeno(x)
    alleles <- unlist(x$data[1, , drop = TRUE])
    alleles <- alleles[!is.na(alleles) & nzchar(alleles)]
    return(unique(as.character(alleles)))
  }
  if (is.character(x)) {
    # Check if these look like alleles (contain *)
    if (all(grepl("\\*", x))) {
      return(unique(x))
    }
  }
  stop("Input must be an hla_genotype object or character vector of allele names.")
}


#' Get peptides from donor input
#' @keywords internal
.getPeptides <- function(donor, recipient, peptide_length) {
  # If donor is already peptides (character without *)
  if (is.character(donor) && !any(grepl("\\*", donor))) {
    # Filter to requested peptide length
    peptides <- donor[nchar(donor) == peptide_length]
    return(unique(peptides))
  }

  # Otherwise, derive peptides from sequence comparison
  donor_alleles <- .extractAlleles(donor)
  recipient_alleles <- .extractAlleles(recipient)

  # Get mismatched positions and generate peptides
  peptides <- .deriveMismatchedPeptides(
    donor_alleles = donor_alleles,
    recipient_alleles = recipient_alleles,
    peptide_length = peptide_length
  )
  return(unique(peptides))
}


#' Derive mismatched peptides from sequence comparison
#' @keywords internal
.deriveMismatchedPeptides <- function(donor_alleles, recipient_alleles, peptide_length) {
  peptides <- character(0)

  for (d_allele in donor_alleles) {
    # Get donor sequence
    d_seq <- tryCatch(
      getAlleleSequence(d_allele),
      error = function(e) NULL
    )
    if (is.null(d_seq) || nchar(d_seq) < peptide_length) next

    # Compare against each recipient allele of same locus
    d_locus <- sub("\\*.*", "", d_allele)

    for (r_allele in recipient_alleles) {
      r_locus <- sub("\\*.*", "", r_allele)
      if (d_locus != r_locus) next

      r_seq <- tryCatch(
        getAlleleSequence(r_allele),
        error = function(e) NULL
      )
      if (is.null(r_seq)) next

      # Find mismatch positions
      mismatch_detail <- quantifyMismatch(r_seq, d_seq, return = "detail")
      # Filter to only actual mismatches
      mismatch_detail <- mismatch_detail[mismatch_detail$is_mismatch, , drop = FALSE]
      if (nrow(mismatch_detail) == 0) next

      # Generate peptides around mismatch positions
      mismatch_positions <- mismatch_detail$alignment_position

      for (pos in mismatch_positions) {
        # Generate all peptides that include this position
        start_min <- max(1, pos - peptide_length + 1)
        start_max <- min(pos, nchar(d_seq) - peptide_length + 1)

        if (start_max >= start_min) {
          for (start in start_min:start_max) {
            pep <- substr(d_seq, start, start + peptide_length - 1)
            if (nchar(pep) == peptide_length && !grepl("-|X", pep)) {
              peptides <- c(peptides, pep)
            }
          }
        }
      }
    }
  }

  return(peptides)
}


#' Predict binding using selected backend
#' @keywords internal
.predictBinding <- function(peptides, alleles, backend, backend_path, peptide_length) {
  if (backend == "pwm") {
    return(.predictBindingPWM(peptides, alleles))
  } else if (backend == "netmhcpan") {
    return(.predictBindingNetMHCpan(peptides, alleles, backend_path))
  } else if (backend == "mhcnuggets") {
    return(.predictBindingMHCnuggets(peptides, alleles))
  }
  stop("Unknown backend: ", backend)
}


#' PWM-based binding prediction (built-in)
#' @keywords internal
.predictBindingPWM <- function(peptides, alleles) {
  # Simplified PWM-based prediction using HLA supertype scoring
  # This provides a reasonable approximation without external dependencies

  # HLA supertype mapping (simplified)
  supertypes <- list(
    A02 = c("A*02:01", "A*02:02", "A*02:03", "A*02:06", "A*02:07", "A*02:11", "A*68:02"),
    A03 = c("A*03:01", "A*11:01", "A*31:01", "A*33:01", "A*68:01"),
    A01 = c("A*01:01", "A*26:01", "A*32:01"),
    A24 = c("A*24:02", "A*23:01"),
    B07 = c("B*07:02", "B*35:01", "B*51:01", "B*53:01", "B*54:01", "B*55:01"),
    B44 = c("B*44:02", "B*44:03", "B*18:01", "B*37:01", "B*40:01", "B*41:01"),
    B27 = c("B*27:05", "B*14:01", "B*38:01", "B*39:01"),
    B58 = c("B*58:01", "B*57:01"),
    B08 = c("B*08:01")
  )

  # Preferred anchor residues for each supertype (position 2 and C-terminus)
  # These are simplified for the PWM approach
  anchor_p2 <- list(
    A02 = c("L", "M", "V", "I", "A", "T"),
    A03 = c("L", "V", "M", "I", "S", "A", "T"),
    A01 = c("T", "S", "M", "L"),
    A24 = c("Y", "F", "W"),
    B07 = c("P"),
    B44 = c("E", "D"),
    B27 = c("R", "K", "H"),
    B58 = c("A", "S", "T"),
    B08 = c("K", "R")
  )

  anchor_c <- list(
    A02 = c("L", "V", "I", "M", "A"),
    A03 = c("K", "R", "Y"),
    A01 = c("Y"),
    A24 = c("F", "L", "I", "W"),
    B07 = c("L", "M", "F"),
    B44 = c("Y", "F", "W"),
    B27 = c("L", "F", "K", "R"),
    B58 = c("W", "F", "Y"),
    B08 = c("L")
  )

  results <- data.frame(
    peptide = character(0),
    hla_allele = character(0),
    predicted_ic50 = numeric(0),
    stringsAsFactors = FALSE
  )

  for (allele in alleles) {
    # Find supertype for this allele
    supertype <- NULL
    for (st in names(supertypes)) {
      # Check exact match or prefix match
      if (allele %in% supertypes[[st]] ||
          any(startsWith(allele, sub("\\*.*", "*", supertypes[[st]])))) {
        supertype <- st
        break
      }
    }

    # Default to A02-like if unknown
    if (is.null(supertype)) supertype <- "A02"

    for (pep in peptides) {
      if (nchar(pep) < 2) next

      # Score based on anchor positions
      p2 <- substr(pep, 2, 2)
      pc <- substr(pep, nchar(pep), nchar(pep))

      score <- 0

      # P2 anchor contribution
      if (p2 %in% anchor_p2[[supertype]]) {
        score <- score + 2
      } else if (p2 %in% c("L", "M", "V", "I", "A")) {
        score <- score + 1
      }

      # C-terminal anchor contribution
      if (pc %in% anchor_c[[supertype]]) {
        score <- score + 2
      } else if (pc %in% c("L", "V", "I", "F", "Y")) {
        score <- score + 1
      }

      # Convert score to approximate IC50
      # Higher score = lower IC50 (better binding)
      # Score 4 -> ~100nM, Score 0 -> ~10000nM
      ic50 <- 10000 / (2^score)
      ic50 <- max(50, min(50000, ic50))

      results <- rbind(results, data.frame(
        peptide = pep,
        hla_allele = allele,
        predicted_ic50 = ic50,
        stringsAsFactors = FALSE
      ))
    }
  }

  return(results)
}


#' NetMHCpan-based binding prediction
#' @keywords internal
.predictBindingNetMHCpan <- function(peptides, alleles, backend_path) {
  if (is.null(backend_path) || !file.exists(backend_path)) {
    stop("NetMHCpan backend requires 'backend_path' to point to the netMHCpan executable.")
  }

  results <- data.frame(
    peptide = character(0),
    hla_allele = character(0),
    predicted_ic50 = numeric(0),
    stringsAsFactors = FALSE
  )

  # Create temp file for peptides
  pep_file <- tempfile(fileext = ".pep")
  writeLines(peptides, pep_file)
  on.exit(unlink(pep_file), add = TRUE)

  for (allele in alleles) {
    # Format allele for NetMHCpan (e.g., HLA-A02:01)
    formatted_allele <- gsub("\\*", "", allele)
    formatted_allele <- paste0("HLA-", formatted_allele)

    # Run NetMHCpan
    out_file <- tempfile(fileext = ".out")
    on.exit(unlink(out_file), add = TRUE)

    cmd <- sprintf(
      "%s -p %s -a %s -BA > %s 2>&1",
      backend_path, pep_file, formatted_allele, out_file
    )

    system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)

    # Parse output
    if (file.exists(out_file)) {
      lines <- readLines(out_file, warn = FALSE)
      # Find data lines (they start with position numbers after header)
      data_lines <- grep("^\\s*[0-9]", lines, value = TRUE)

      for (line in data_lines) {
        fields <- strsplit(trimws(line), "\\s+")[[1]]
        if (length(fields) >= 12) {
          pep <- fields[3]
          ic50 <- as.numeric(fields[12])
          if (!is.na(ic50) && pep %in% peptides) {
            results <- rbind(results, data.frame(
              peptide = pep,
              hla_allele = allele,
              predicted_ic50 = ic50,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }

  # If NetMHCpan didn't return results for some peptides, add them with NA
  missing <- setdiff(
    paste(rep(peptides, each = length(alleles)), rep(alleles, length(peptides)), sep = "_"),
    paste(results$peptide, results$hla_allele, sep = "_")
  )

  if (length(missing) > 0) {
    missing_parts <- strsplit(missing, "_")
    for (mp in missing_parts) {
      results <- rbind(results, data.frame(
        peptide = mp[1],
        hla_allele = mp[2],
        predicted_ic50 = 50000,  # Assume non-binder
        stringsAsFactors = FALSE
      ))
    }
  }

  return(results)
}


#' MHCnuggets-based binding prediction
#'
#' Uses the predictMHCnuggets function for deep learning-based binding prediction.
#' MHCnuggets is available at \url{https://github.com/KarchinLab/mhcnuggets}.
#'
#' @keywords internal
.predictBindingMHCnuggets <- function(peptides, alleles) {
  results <- data.frame(
    peptide = character(0),
    hla_allele = character(0),
    predicted_ic50 = numeric(0),
    stringsAsFactors = FALSE
  )


  # Run predictions for each allele
  for (allele in alleles) {
    # Determine MHC class from allele name
    locus <- sub("\\*.*", "", allele)
    mhc_class <- if (locus %in% c("A", "B", "C")) "I" else "II"

    tryCatch({
      # Use the package's predictMHCnuggets function
      pred_result <- predictMHCnuggets(
        peptides = peptides,
        allele = allele,
        mhc_class = mhc_class
      )

      results <- rbind(results, data.frame(
        peptide = pred_result$peptide,
        hla_allele = allele,
        predicted_ic50 = pred_result$ic50,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      # If prediction fails for this allele, add entries with high IC50
      warning(sprintf("MHCnuggets prediction failed for allele %s: %s", allele, e$message))
      results <<- rbind(results, data.frame(
        peptide = peptides,
        hla_allele = allele,
        predicted_ic50 = 50000,
        stringsAsFactors = FALSE
      ))
    })
  }

  return(results)
}


#' Visualize Cross-Locus Peptide Binding Results
#'
#' @description
#' Creates visualizations of peptide binding predictions across all loci.
#' This function is designed for advanced cross-locus analysis where peptides
#' from multiple donor alleles are tested against multiple recipient alleles.
#'
#' @param binding_results A list containing an `all_predictions` data.frame with columns:
#'   `donor_allele`, `recipient_allele`, `binding` (logical), `recipient_locus`,
#'   `mhc_class`, `donor_locus`, and optionally `ic50`.
#' @param plot_type Type of plot: "heatmap", "bar_by_recipient", "bar_by_donor", or "scatter"
#' @param palette Character. A color palette name. Defaults to "spectral".
#' @param ... Additional arguments passed to the ggplot theme.
#'
#' @return ggplot object
#'
#' @examples
#' # Create example binding results data structure
#' binding_results <- list(
#'   all_predictions = data.frame(
#'     donor_allele = rep(c("A*01:01", "A*24:02"), each = 4),
#'     recipient_allele = rep(c("A*02:01", "A*03:01"), 4),
#'     recipient_locus = "A",
#'     donor_locus = "A",
#'     mhc_class = "I",
#'     binding = c(TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
#'     ic50 = c(100, 6000, 250, 150, 8000, 300, 7500, 9000)
#'   )
#' )
#'
#' # Create heatmap visualization
#' p <- visualizePeptideBinding(binding_results, plot_type = "heatmap")
#' print(p)
#'
#' # Create bar plot by recipient
#' p2 <- visualizePeptideBinding(binding_results, plot_type = "bar_by_recipient")
#' print(p2)
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

# Global variable declarations to avoid R CMD check notes
utils::globalVariables(c(
 "xmin", "xmax", "exon_num", "midpoint", "label",
 "Sequence", "Residue", "Position",
 "x_min", "x_max", "y_min", "y_max",
 "PosLabel", "x", "color", "aas", "group"
))

#' Plot HLA Sequence Alignments with Enhanced Visualization
#'
#' @description
#' Compares two or more HLA allele sequences (nucleotide or protein) with
#' comprehensive visualization including domain/exon boundaries, mismatch
#' highlighting, and biophysical property annotations. Uses `getAlleleSequence()`
#' for sequence retrieval and `quantifyMismatch()` for mismatch analysis.
#'
#' @param alleles Character vector of HLA allele names to compare (2 or more).
#'   Names should be in standard format (e.g., "A*01:01", "DRB1*03:01").
#' @param sequences Optional named list or character vector of sequences. If provided,
#'   these sequences are used instead of fetching from IMGT/HLA. Names should
#'   correspond to `alleles` or be descriptive labels.
#' @param locus HLA locus for exon boundary annotation (e.g., "A", "B", "DRB1").
#'   Only needed when `show_exons = TRUE`. If NULL, will be inferred from allele names.
#' @param seq_type Type of sequence: "protein" or "nucleotide" (default "protein").
#' @param focus_on_mismatches If TRUE, show zoomed view of mismatch regions (default FALSE).
#' @param context_window Number of positions to show around mismatches when
#'   `focus_on_mismatches = TRUE` (default 20).
#' @param show_exons If TRUE, display exon/domain boundaries (default TRUE).
#' @param color_scheme For protein sequences: "properties" (physicochemical grouping),
#'   "hydropathy" (hydrophobicity scale), or "classic" (rainbow). Default "properties".
#' @param reference_idx Integer index of the sequence to use as reference for
#'   mismatch counting (default 1, the first sequence).
#' @param filter_charge Passed to `quantifyMismatch()`. NULL (default), TRUE, or FALSE.
#' @param filter_polarity Passed to `quantifyMismatch()`. NULL (default), TRUE, or FALSE.
#' @param use_cache Whether to use caching for sequence retrieval (default TRUE).
#' @param verbose Print progress messages (default TRUE).
#'
#' @return A list with class "hla_sequence_plot" containing:
#' \describe{
#'   \item{plot}{Combined ggplot2 plot object using patchwork}
#'   \item{mismatch_summary}{Data frame summarizing pairwise mismatches}
#'   \item{mismatch_details}{List of detailed mismatch information for each comparison}
#'   \item{exon_boundaries}{Numeric vector of exon/domain boundaries (if available)}
#'   \item{exon_descriptions}{Character vector describing each domain/exon}
#'   \item{alignment_df}{Data frame of aligned sequences for custom plotting}
#'   \item{sequences}{Named list of the sequences used}
#'   \item{plots}{List of individual plot components}
#' }
#'
#' @examples
#' \donttest{
#' # Compare two protein sequences
#' result <- plotHLASequences(
#'   alleles = c("A*01:01", "A*02:01"),
#'   seq_type = "protein"
#' )
#' print(result)
#'
#' # Compare three alleles with charge filtering
#' result <- plotHLASequences(
#'   alleles = c("DRB1*03:01", "DRB1*04:01", "DRB1*07:01"),
#'   seq_type = "protein",
#'   filter_charge = TRUE
#' )
#'
#' # Use custom sequences
#' seqs <- list(
#'   "Seq1" = "YFAMYGEKVAHTHVDTLYVRYHY",
#'   "Seq2" = "YFDMYGEKVAHTHVDTLYVRFHY",
#'   "Seq3" = "YFAMYGEKVAHTHVDTLYVRFHY"
#' )
#' result <- plotHLASequences(
#'   alleles = names(seqs),
#'   sequences = seqs,
#'   seq_type = "protein",
#'   show_exons = FALSE
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_rect geom_text geom_vline
#'   geom_segment scale_fill_manual scale_color_manual scale_x_continuous
#'   labs theme_minimal theme element_text element_blank margin coord_cartesian
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom stats setNames
#' @export
plotHLASequences <- function(alleles,
                              sequences = NULL,
                              locus = NULL,
                              seq_type = c("protein", "nucleotide"),
                              focus_on_mismatches = FALSE,
                              context_window = 20,
                              show_exons = TRUE,
                              color_scheme = c("properties", "hydropathy", "classic"),
                              reference_idx = 1,
                              filter_charge = NULL,
                              filter_polarity = NULL,
                              use_cache = TRUE,
                              verbose = TRUE) {

  # Validate arguments
  seq_type <- match.arg(seq_type)
  color_scheme <- match.arg(color_scheme)

  if (length(alleles) < 2) {
    stop("At least 2 alleles/sequences are required for comparison.")
  }

  if (!is.null(reference_idx) && (reference_idx < 1 || reference_idx > length(alleles))) {
    stop("reference_idx must be between 1 and the number of alleles.")
  }

  # ==================================================================
  # 1. RETRIEVE OR VALIDATE SEQUENCES
  # ==================================================================

  if (is.null(sequences)) {
    # Fetch sequences using getAlleleSequence
    type_arg <- ifelse(seq_type == "protein", "PROT", "NUC")

    if (verbose) message("Retrieving sequences from IMGT/HLA...")

    sequences <- lapply(alleles, function(a) {
      tryCatch(
        getAlleleSequence(a, type = type_arg, use_cache = use_cache),
        error = function(e) {
          stop(sprintf("Failed to retrieve sequence for '%s': %s", a, e$message))
        }
      )
    })
    names(sequences) <- alleles

    if (verbose) message(sprintf("Retrieved %d sequences.", length(sequences)))
  } else {
    # Validate provided sequences
    if (is.character(sequences) && !is.null(names(sequences))) {
      sequences <- as.list(sequences)
    }

    if (!is.list(sequences)) {
      stop("sequences must be a named list or named character vector.")
    }

    if (length(sequences) != length(alleles)) {
      stop("Length of sequences must match length of alleles.")
    }

    # Use alleles as names if sequences aren't named
    if (is.null(names(sequences))) {
      names(sequences) <- alleles
    }
  }

  # Infer locus from allele names if not provided
  if (is.null(locus) && show_exons) {
    locus <- .inferLocus(alleles[1])
    if (verbose && !is.null(locus)) {
      message(sprintf("Inferred locus: %s", locus))
    }
  }

  # ==================================================================
  # 2. PERFORM MULTIPLE SEQUENCE ALIGNMENT
  # ==================================================================

  alignment_result <- .performMSA(sequences, seq_type, verbose)
  aligned_seqs <- alignment_result$aligned
  alignment_length <- alignment_result$length

  # ==================================================================
  # 3. GET EXON BOUNDARIES
  # ==================================================================

  exon_info <- NULL
  exon_boundaries <- NULL
  exon_descriptions <- NULL

  if (show_exons && !is.null(locus)) {
    exon_info <- .getKnownExonBoundaries(locus, seq_type)
    if (!is.null(exon_info)) {
      exon_boundaries <- exon_info$boundaries
      exon_descriptions <- exon_info$description

      # Filter boundaries within alignment
      valid_idx <- exon_boundaries < alignment_length
      exon_boundaries <- exon_boundaries[valid_idx]

      if (verbose && length(exon_boundaries) > 0) {
        message(sprintf("Using exon boundaries for %s: %s",
                        locus, paste(exon_boundaries, collapse = ", ")))
      }
    }
  }

  # ==================================================================
  # 4. COMPUTE PAIRWISE MISMATCHES
  # ==================================================================

  n_seqs <- length(aligned_seqs)
  seq_names <- names(aligned_seqs)
  ref_name <- seq_names[reference_idx]
  ref_seq <- aligned_seqs[[reference_idx]]

  mismatch_summary <- data.frame(
    comparison = character(),
    total_mismatches = integer(),
    filtered_mismatches = integer(),
    comparable_positions = integer(),
    percent_identity = numeric(),
    stringsAsFactors = FALSE
  )

  mismatch_details <- list()
  all_mismatch_positions <- c()

  for (i in seq_len(n_seqs)) {
    if (i == reference_idx) next

    comp_name <- seq_names[i]
    comp_seq <- aligned_seqs[[i]]

    if (seq_type == "protein") {
      # Use quantifyMismatch for protein sequences
      # Get total mismatches
      total_mm <- quantifyMismatch(
        ref_seq, comp_seq,
        filter_charge = NULL,
        filter_polarity = NULL,
        return = "count",
        na_action = "exclude",
        count_gaps = TRUE
      )

      # Get filtered mismatches
      filtered_mm <- quantifyMismatch(
        ref_seq, comp_seq,
        filter_charge = filter_charge,
        filter_polarity = filter_polarity,
        return = "count",
        na_action = "exclude",
        count_gaps = TRUE
      )

      # Get detailed breakdown
      details <- quantifyMismatch(
        ref_seq, comp_seq,
        filter_charge = filter_charge,
        filter_polarity = filter_polarity,
        return = "detail",
        na_action = "exclude",
        count_gaps = TRUE
      )

      # Calculate comparable positions and identity
      non_gap_positions <- sum(!details$is_gap_ref & !details$is_gap_alt)
      pct_identity <- round(100 * (non_gap_positions - total_mm) / non_gap_positions, 2)

      # Record mismatch positions
      mm_pos <- details$alignment_position[details$is_mismatch]

      # Add domain info if available
      if (!is.null(exon_boundaries)) {
        details$domain <- sapply(details$alignment_position, function(pos) {
          .getDomain(pos, exon_boundaries, exon_descriptions)
        })
      }
    } else {
      # For nucleotide sequences, use simple character comparison
      details <- .compareNucleotideSequences(ref_seq, comp_seq)

      total_mm <- sum(details$is_mismatch)
      filtered_mm <- total_mm  # No filtering for nucleotides
      non_gap_positions <- sum(!details$is_gap_ref & !details$is_gap_alt)
      pct_identity <- round(100 * (non_gap_positions - total_mm) / max(non_gap_positions, 1), 2)

      # Record mismatch positions
      mm_pos <- details$alignment_position[details$is_mismatch]

      # Add domain info if available
      if (!is.null(exon_boundaries)) {
        details$domain <- sapply(details$alignment_position, function(pos) {
          .getDomain(pos, exon_boundaries, exon_descriptions)
        })
      }
    }

    all_mismatch_positions <- c(all_mismatch_positions, mm_pos)
    mismatch_details[[paste(ref_name, "vs", comp_name)]] <- details

    mismatch_summary <- rbind(mismatch_summary, data.frame(
      comparison = paste(ref_name, "vs", comp_name),
      total_mismatches = total_mm,
      filtered_mismatches = filtered_mm,
      comparable_positions = non_gap_positions,
      percent_identity = pct_identity,
      stringsAsFactors = FALSE
    ))
  }

  all_mismatch_positions <- sort(unique(all_mismatch_positions))

  if (verbose) {
    message(sprintf("Found %d unique mismatch positions across all comparisons.",
                    length(all_mismatch_positions)))
  }

  # ==================================================================
  # 5. CREATE DATA FRAMES FOR PLOTTING
  # ==================================================================

  # Split aligned sequences into characters
  seq_chars <- lapply(aligned_seqs, function(s) strsplit(s, "")[[1]])

  # Create alignment data frame
  align_df <- do.call(rbind, lapply(seq_along(seq_chars), function(i) {
    data.frame(
      Position = seq_len(alignment_length),
      Sequence = seq_names[i],
      Residue = seq_chars[[i]],
      stringsAsFactors = FALSE
    )
  }))

  # Set sequence factor levels (reverse for plotting)
  align_df$Sequence <- factor(align_df$Sequence, levels = rev(seq_names))

  # Get color palette
  pal <- .getColorPalette(seq_type, color_scheme)

  # Create exon data frame
  exon_df <- NULL
  if (!is.null(exon_boundaries) && length(exon_boundaries) > 0) {
    exon_starts <- c(1, exon_boundaries + 1)
    exon_ends <- c(exon_boundaries, alignment_length)

    exon_df <- data.frame(
      xmin = exon_starts - 0.5,
      xmax = exon_ends + 0.5,
      exon_num = seq_along(exon_starts),
      label = exon_descriptions[seq_along(exon_starts)],
      stringsAsFactors = FALSE
    )
    exon_df$midpoint <- (exon_df$xmin + exon_df$xmax) / 2
  }

  # Mismatch rectangles
  mismatch_rect <- NULL
  if (length(all_mismatch_positions) > 0) {
    mismatch_rect <- data.frame(
      x_min = all_mismatch_positions - 0.5,
      x_max = all_mismatch_positions + 0.5,
      y_min = 0.5,
      y_max = n_seqs + 0.5
    )
  }

  # ==================================================================
  # 6. CREATE PLOTS
  # ==================================================================

  plots <- list()

  # === PLOT A: Exon/Domain Map Bar ===
  if (!is.null(exon_df)) {
    plots$exon_bar <- ggplot2::ggplot(exon_df) +
      ggplot2::geom_rect(
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1,
                     fill = factor(exon_num %% 2)),
        color = "black", linewidth = 0.4
      ) +
      ggplot2::geom_text(
        ggplot2::aes(x = midpoint, y = 0.5, label = label),
        size = 3, fontface = "bold"
      ) +
      ggplot2::scale_fill_manual(
        values = c("0" = "#e2e8f0", "1" = "#f7fafc"), guide = "none"
      ) +
      ggplot2::labs(title = "Domain Structure") +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 11, hjust = 0),
        plot.margin = ggplot2::margin(5, 10, 5, 10)
      ) +
      ggplot2::coord_cartesian(xlim = c(0.5, alignment_length + 0.5))
  }

  # === PLOT B: Main Alignment ===
  plots$main <- ggplot2::ggplot(align_df, ggplot2::aes(x = Position, y = Sequence)) +
    {if (!is.null(exon_boundaries)) {
      ggplot2::geom_vline(
        xintercept = exon_boundaries + 0.5,
        linetype = "dashed", color = "grey40", linewidth = 0.5
      )
    }} +
    ggplot2::geom_tile(ggplot2::aes(fill = Residue), height = 0.9) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    {if (!is.null(mismatch_rect) && nrow(mismatch_rect) > 0) {
      ggplot2::geom_rect(
        data = mismatch_rect, inherit.aes = FALSE,
        ggplot2::aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max),
        color = "red", fill = NA, linewidth = 1.2
      )
    }} +
    ggplot2::labs(x = "Alignment Position", y = NULL, title = "Full Alignment") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold", size = 9),
      plot.title = ggplot2::element_text(face = "bold", size = 11)
    )

  # === PLOT C: Zoomed View (if mismatches exist) ===
  if (length(all_mismatch_positions) > 0 && focus_on_mismatches) {
    zoom_start <- max(1, min(all_mismatch_positions) - context_window)
    zoom_end <- min(alignment_length, max(all_mismatch_positions) + context_window)

    zoom_df <- align_df[align_df$Position >= zoom_start &
                          align_df$Position <= zoom_end, ]
    zoom_mismatch_rect <- NULL
    if (!is.null(mismatch_rect)) {
      zoom_mismatch_rect <- mismatch_rect[
        mismatch_rect$x_min >= (zoom_start - 0.5) &
          mismatch_rect$x_max <= (zoom_end + 0.5),
      ]
    }

    zoom_exon_bounds <- NULL
    if (!is.null(exon_boundaries)) {
      zoom_exon_bounds <- exon_boundaries[
        exon_boundaries >= zoom_start & exon_boundaries <= zoom_end
      ]
    }

    plots$zoom <- ggplot2::ggplot(zoom_df, ggplot2::aes(x = Position, y = Sequence)) +
      {if (length(zoom_exon_bounds) > 0) {
        ggplot2::geom_vline(
          xintercept = zoom_exon_bounds + 0.5,
          linetype = "dashed", color = "grey40", linewidth = 0.5
        )
      }} +
      ggplot2::geom_tile(
        ggplot2::aes(fill = Residue),
        height = 0.85, color = "white", linewidth = 0.2
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = Residue),
        size = 2.8, color = "white", fontface = "bold"
      ) +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      {if (!is.null(zoom_mismatch_rect) && nrow(zoom_mismatch_rect) > 0) {
        ggplot2::geom_rect(
          data = zoom_mismatch_rect, inherit.aes = FALSE,
          ggplot2::aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max),
          color = "red", fill = NA, linewidth = 1.5
        )
      }} +
      ggplot2::labs(
        x = "Position", y = NULL,
        title = sprintf("Zoomed View (positions %d-%d)", zoom_start, zoom_end)
      ) +
      ggplot2::scale_x_continuous(breaks = seq(zoom_start, zoom_end, by = 5)) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(face = "bold", size = 9),
        axis.text.x = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(face = "bold", size = 11)
      )
  }

  # === PLOT D: Mismatch Details Grid ===
  if (length(all_mismatch_positions) > 0) {
    mismatch_zoom_df <- do.call(rbind, lapply(seq_along(seq_chars), function(i) {
      data.frame(
        Position = all_mismatch_positions,
        Sequence = seq_names[i],
        Residue = seq_chars[[i]][all_mismatch_positions],
        stringsAsFactors = FALSE
      )
    }))
    mismatch_zoom_df$Sequence <- factor(mismatch_zoom_df$Sequence, levels = rev(seq_names))
    mismatch_zoom_df$PosLabel <- factor(
      mismatch_zoom_df$Position,
      levels = sort(unique(all_mismatch_positions))
    )

    plots$mismatch <- ggplot2::ggplot(
      mismatch_zoom_df,
      ggplot2::aes(x = PosLabel, y = Sequence)
    ) +
      ggplot2::geom_tile(
        ggplot2::aes(fill = Residue),
        height = 0.85, color = "black", linewidth = 0.5
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = Residue),
        size = 4, fontface = "bold", color = "white"
      ) +
      ggplot2::scale_fill_manual(values = pal, guide = "none") +
      ggplot2::labs(x = "Position", y = NULL, title = "Mismatch Positions") +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(face = "bold", size = 10),
        axis.text.x = ggplot2::element_text(size = 10, face = "bold"),
        plot.title = ggplot2::element_text(face = "bold", size = 11)
      )
  }

  # === PLOT E: Property Legend (protein only) ===
  if (seq_type == "protein" && color_scheme == "properties") {
    prop_groups <- .getPropertyGroups()
    legend_df <- data.frame(
      group = names(prop_groups),
      color = sapply(prop_groups, function(x) x$color),
      aas = sapply(prop_groups, function(x) paste(x$aas, collapse = ", ")),
      x = seq_along(prop_groups),
      stringsAsFactors = FALSE
    )

    plots$legend <- ggplot2::ggplot(legend_df, ggplot2::aes(x = x, y = 1)) +
      ggplot2::geom_tile(
        ggplot2::aes(fill = I(color)),
        height = 0.7, width = 0.9
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = aas),
        y = 1, size = 2.5, color = "white", fontface = "bold"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = group),
        y = 0.35, size = 2.2, lineheight = 0.9
      ) +
      ggplot2::labs(title = "Amino Acid Properties") +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 10, hjust = 0.5),
        plot.margin = ggplot2::margin(5, 5, 10, 5)
      ) +
      ggplot2::coord_cartesian(ylim = c(-0.1, 1.5), clip = "off")
  }

  # ==================================================================
  # 7. COMPOSE FINAL PLOT
  # ==================================================================

  # Build title and subtitle
  title_text <- sprintf("HLA Sequence Comparison (%d sequences)", n_seqs)

  total_mm <- sum(mismatch_summary$filtered_mismatches)
  subtitle_parts <- sprintf(
    "%d total mismatch%s (reference: %s)",
    total_mm,
    ifelse(total_mm != 1, "es", ""),
    ref_name
  )

  if (!is.null(filter_charge) || !is.null(filter_polarity)) {
    filters <- c()
    if (!is.null(filter_charge)) {
      filters <- c(filters, sprintf("charge=%s", filter_charge))
    }
    if (!is.null(filter_polarity)) {
      filters <- c(filters, sprintf("polarity=%s", filter_polarity))
    }
    subtitle_parts <- paste0(subtitle_parts, " | Filters: ", paste(filters, collapse = ", "))
  }

  # Combine plots with patchwork
  plots_to_combine <- list()
  heights <- c()

  if (!is.null(plots$exon_bar)) {
    plots_to_combine <- c(plots_to_combine, list(plots$exon_bar))
    heights <- c(heights, 0.8)
  }

  plots_to_combine <- c(plots_to_combine, list(plots$main))
  heights <- c(heights, max(2, 0.8 * n_seqs))

  if (!is.null(plots$zoom)) {
    plots_to_combine <- c(plots_to_combine, list(plots$zoom))
    heights <- c(heights, max(1.5, 0.6 * n_seqs))
  }

  if (!is.null(plots$mismatch)) {
    plots_to_combine <- c(plots_to_combine, list(plots$mismatch))
    heights <- c(heights, max(1.5, 0.5 * n_seqs))
  }

  if (!is.null(plots$legend)) {
    plots_to_combine <- c(plots_to_combine, list(plots$legend))
    heights <- c(heights, 1.2)
  }

  if (length(plots_to_combine) > 1) {
    final_plot <- patchwork::wrap_plots(plots_to_combine, ncol = 1, heights = heights) +
      patchwork::plot_annotation(
        title = title_text,
        subtitle = subtitle_parts,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 14),
          plot.subtitle = ggplot2::element_text(size = 11, color = "grey40")
        )
      )
  } else {
    final_plot <- plots$main +
      ggplot2::labs(title = title_text, subtitle = subtitle_parts)
  }

  # ==================================================================
  # 8. RETURN RESULTS
  # ==================================================================

  result <- list(
    plot = final_plot,
    mismatch_summary = mismatch_summary,
    mismatch_details = mismatch_details,
    exon_boundaries = exon_boundaries,
    exon_descriptions = exon_descriptions,
    alignment_df = align_df,
    sequences = sequences,
    aligned_sequences = aligned_seqs,
    plots = plots
  )

  class(result) <- c("hla_sequence_plot", "list")
  return(result)
}

#' @export
print.hla_sequence_plot <- function(x, ...) {
  print(x$plot)
  invisible(x)
}

# ==================================================================
# INTERNAL HELPER FUNCTIONS
# ==================================================================

#' Infer locus from allele name
#' @noRd
.inferLocus <- function(allele) {
  # Extract locus from allele name like "A*01:01" or "DRB1*03:01"
  match <- regmatches(allele, regexpr("^[A-Z]+[0-9]*", allele))
  if (length(match) > 0 && nchar(match) > 0) {
    return(match)
  }
  return(NULL)
}

#' Get domain name for a position
#' @noRd
.getDomain <- function(pos, boundaries, descriptions) {
  for (i in seq_along(boundaries)) {
    if (pos <= boundaries[i]) {
      return(descriptions[i])
    }
  }
  return(descriptions[length(descriptions)])
}

#' Perform multiple sequence alignment
#' @noRd
.performMSA <- function(sequences, seq_type, verbose = FALSE) {
  n_seqs <- length(sequences)

  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings package is required for sequence alignment.")
  }

  if (!requireNamespace("pwalign", quietly = TRUE)) {
    stop("pwalign package is required for sequence alignment.")
  }

  # For 2 sequences, use pairwise alignment
  if (n_seqs == 2) {
    if (seq_type == "protein") {
      data(list = "BLOSUM62", package = "pwalign", envir = environment())
      mat <- get("BLOSUM62")

      pwa <- pwalign::pairwiseAlignment(
        pattern = Biostrings::AAString(sequences[[1]]),
        subject = Biostrings::AAString(sequences[[2]]),
        type = "global",
        substitutionMatrix = mat,
        gapOpening = 10,
        gapExtension = 0.5
      )
    } else {
      # Nucleotide alignment
      dna_letters <- c("A", "C", "G", "T")
      custom_mat <- matrix(0, nrow = 4, ncol = 4,
                           dimnames = list(dna_letters, dna_letters))
      diag(custom_mat) <- 1

      pwa <- pwalign::pairwiseAlignment(
        pattern = Biostrings::DNAString(sequences[[1]]),
        subject = Biostrings::DNAString(sequences[[2]]),
        type = "global",
        substitutionMatrix = custom_mat,
        gapOpening = 10,
        gapExtension = 4
      )
    }

    aligned <- list(
      as.character(pwalign::pattern(pwa)),
      as.character(pwalign::subject(pwa))
    )
    names(aligned) <- names(sequences)

    return(list(
      aligned = aligned,
      length = nchar(aligned[[1]])
    ))
  }

  # For 3+ sequences, perform progressive alignment using first as reference
  if (verbose) message("Performing progressive alignment with reference sequence...")

  ref_seq <- sequences[[1]]
  aligned <- list()
  aligned[[names(sequences)[1]]] <- ref_seq

  # Track alignment positions
  ref_aligned <- ref_seq

  for (i in 2:n_seqs) {
    if (seq_type == "protein") {
      data(list = "BLOSUM62", package = "pwalign", envir = environment())
      mat <- get("BLOSUM62")

      pwa <- pwalign::pairwiseAlignment(
        pattern = Biostrings::AAString(ref_aligned),
        subject = Biostrings::AAString(sequences[[i]]),
        type = "global",
        substitutionMatrix = mat,
        gapOpening = 10,
        gapExtension = 0.5
      )
    } else {
      dna_letters <- c("A", "C", "G", "T")
      custom_mat <- matrix(0, nrow = 4, ncol = 4,
                           dimnames = list(dna_letters, dna_letters))
      diag(custom_mat) <- 1

      pwa <- pwalign::pairwiseAlignment(
        pattern = Biostrings::DNAString(ref_aligned),
        subject = Biostrings::DNAString(sequences[[i]]),
        type = "global",
        substitutionMatrix = custom_mat,
        gapOpening = 10,
        gapExtension = 4
      )
    }

    new_ref <- as.character(pwalign::pattern(pwa))
    new_sub <- as.character(pwalign::subject(pwa))

    # If reference gained gaps, update all previous alignments
    if (nchar(new_ref) > nchar(ref_aligned)) {
      gap_positions <- which(strsplit(new_ref, "")[[1]] == "-")

      for (j in seq_along(aligned)) {
        old_chars <- strsplit(aligned[[j]], "")[[1]]
        new_chars <- character(nchar(new_ref))
        old_idx <- 1

        for (k in seq_len(nchar(new_ref))) {
          if (k %in% gap_positions) {
            new_chars[k] <- "-"
          } else {
            new_chars[k] <- old_chars[old_idx]
            old_idx <- old_idx + 1
          }
        }
        aligned[[j]] <- paste(new_chars, collapse = "")
      }
    }

    ref_aligned <- new_ref
    aligned[[names(sequences)[i]]] <- new_sub
  }

  # Ensure all sequences have same length
  max_len <- max(sapply(aligned, nchar))
  aligned <- lapply(aligned, function(s) {
    if (nchar(s) < max_len) {
      paste0(s, paste(rep("-", max_len - nchar(s)), collapse = ""))
    } else {
      s
    }
  })

  return(list(
    aligned = aligned,
    length = max_len
  ))
}

#' Get known exon boundaries for HLA loci
#' @noRd
.getKnownExonBoundaries <- function(locus, seq_type) {
  # Protein exon boundaries (for mature protein)
  protein_boundaries <- list(
    # Class I molecules
    "A" = list(
      description = c("Leader", "alpha1 domain", "alpha2 domain", "alpha3 domain", "TM+Cyto"),
      boundaries = c(90, 182, 274)
    ),
    "B" = list(
      description = c("Leader", "alpha1 domain", "alpha2 domain", "alpha3 domain", "TM+Cyto"),
      boundaries = c(90, 182, 274)
    ),
    "C" = list(
      description = c("Leader", "alpha1 domain", "alpha2 domain", "alpha3 domain", "TM+Cyto"),
      boundaries = c(90, 182, 274)
    ),

    # Class II alpha chains
    "DRA" = list(
      description = c("Leader", "alpha1 domain", "alpha2 domain", "TM+Cyto"),
      boundaries = c(84, 179)
    ),
    "DQA1" = list(
      description = c("Leader", "alpha1 domain", "alpha2 domain", "TM+Cyto"),
      boundaries = c(87, 181)
    ),
    "DPA1" = list(
      description = c("Leader", "alpha1 domain", "alpha2 domain", "TM+Cyto"),
      boundaries = c(84, 178)
    ),

    # Class II beta chains
    "DRB1" = list(
      description = c("Leader", "beta1 domain", "beta2 domain", "TM+Cyto"),
      boundaries = c(94, 190)
    ),
    "DRB3" = list(
      description = c("Leader", "beta1 domain", "beta2 domain", "TM+Cyto"),
      boundaries = c(94, 190)
    ),
    "DRB4" = list(
      description = c("Leader", "beta1 domain", "beta2 domain", "TM+Cyto"),
      boundaries = c(94, 190)
    ),
    "DRB5" = list(
      description = c("Leader", "beta1 domain", "beta2 domain", "TM+Cyto"),
      boundaries = c(94, 190)
    ),
    "DQB1" = list(
      description = c("Leader", "beta1 domain", "beta2 domain", "TM+Cyto"),
      boundaries = c(94, 188)
    ),
    "DPB1" = list(
      description = c("Leader", "beta1 domain", "beta2 domain", "TM+Cyto"),
      boundaries = c(94, 186)
    )
  )

  # Nucleotide boundaries (CDS positions)
  nucleotide_boundaries <- list(
    "A" = list(
      description = c("Exon 2", "Exon 3", "Exon 4", "Exons 5-8"),
      boundaries = c(270, 546, 822)
    ),
    "B" = list(
      description = c("Exon 2", "Exon 3", "Exon 4", "Exons 5-8"),
      boundaries = c(270, 546, 822)
    ),
    "C" = list(
      description = c("Exon 2", "Exon 3", "Exon 4", "Exons 5-8"),
      boundaries = c(270, 546, 822)
    ),
    "DRB1" = list(
      description = c("Exon 2", "Exon 3", "Exons 4-6"),
      boundaries = c(282, 570)
    ),
    "DQB1" = list(
      description = c("Exon 2", "Exon 3", "Exons 4-6"),
      boundaries = c(282, 564)
    ),
    "DPB1" = list(
      description = c("Exon 2", "Exon 3", "Exons 4-6"),
      boundaries = c(282, 558)
    )
  )

  if (seq_type == "protein") {
    if (locus %in% names(protein_boundaries)) {
      return(protein_boundaries[[locus]])
    }
  } else {
    if (locus %in% names(nucleotide_boundaries)) {
      return(nucleotide_boundaries[[locus]])
    }
  }

  return(NULL)
}

#' Get color palette for sequence visualization
#' @noRd
.getColorPalette <- function(seq_type, color_scheme = "properties") {
  if (seq_type == "nucleotide") {
    return(c(
      "A" = "#E41A1C", "C" = "#377EB8", "G" = "#4DAF4A", "T" = "#FF7F00",
      "U" = "#FF7F00", "N" = "#999999", "-" = "#BDBDBD"
    ))
  }

  if (color_scheme == "properties") {
    return(c(
      # Hydrophobic aliphatic
      "A" = "#2E7D32", "V" = "#388E3C", "I" = "#43A047", "L" = "#4CAF50", "M" = "#66BB6A",
      # Hydrophobic aromatic
      "F" = "#00695C", "W" = "#00796B", "Y" = "#00897B",
      # Polar uncharged
      "S" = "#1565C0", "T" = "#1976D2", "N" = "#1E88E5", "Q" = "#2196F3",
      # Positive charged
      "K" = "#C62828", "R" = "#D32F2F", "H" = "#E53935",
      # Negative charged
      "D" = "#6A1B9A", "E" = "#7B1FA2",
      # Special
      "C" = "#EF6C00", "G" = "#F57C00", "P" = "#FB8C00",
      # Stop and gap
      "*" = "#212121", "-" = "#BDBDBD", "X" = "#9E9E9E"
    ))
  }

  if (color_scheme == "hydropathy") {
    return(c(
      # Most hydrophobic (red)
      "I" = "#B71C1C", "V" = "#C62828", "L" = "#D32F2F",
      "F" = "#E53935", "C" = "#EF5350", "M" = "#EF5350",
      "A" = "#FF7043", "W" = "#FF8A65",
      # Neutral (yellow)
      "G" = "#FDD835", "T" = "#FFEE58", "S" = "#FFF176",
      # Hydrophilic (blue)
      "Y" = "#81D4FA", "P" = "#4FC3F7", "H" = "#29B6F6",
      "N" = "#039BE5", "Q" = "#0288D1", "D" = "#0277BD",
      "E" = "#01579B", "K" = "#1A237E", "R" = "#311B92",
      # Stop and gap
      "*" = "#212121", "-" = "#BDBDBD", "X" = "#9E9E9E"
    ))
  }

  # Classic (rainbow) scheme
  aa_order <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*", "-", "X")
  classic_colors <- c(
    "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
    "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
    "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4", "#b2df8a",
    "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
    "#000000", "#969696", "#9E9E9E"
  )
  return(setNames(classic_colors, aa_order))
}

#' Get amino acid property groups for legend
#' @noRd
.getPropertyGroups <- function() {
  list(
    "Hydrophobic\n(aliphatic)" = list(aas = c("A", "V", "I", "L", "M"), color = "#4CAF50"),
    "Hydrophobic\n(aromatic)" = list(aas = c("F", "W", "Y"), color = "#00796B"),
    "Polar\n(uncharged)" = list(aas = c("S", "T", "N", "Q"), color = "#1976D2"),
    "Positive\n(charged)" = list(aas = c("K", "R", "H"), color = "#D32F2F"),
    "Negative\n(charged)" = list(aas = c("D", "E"), color = "#7B1FA2"),
    "Special" = list(aas = c("C", "G", "P"), color = "#F57C00")
  )
}

#' Compare nucleotide sequences (simple character-by-character)
#' @noRd
.compareNucleotideSequences <- function(seq1, seq2) {
  chars1 <- strsplit(toupper(seq1), "")[[1]]
  chars2 <- strsplit(toupper(seq2), "")[[1]]

  # Ensure same length
  len <- max(length(chars1), length(chars2))
  if (length(chars1) < len) {
    chars1 <- c(chars1, rep("-", len - length(chars1)))
  }
  if (length(chars2) < len) {
    chars2 <- c(chars2, rep("-", len - length(chars2)))
  }

  # Create comparison data frame
  data.frame(
    alignment_position = seq_len(len),
    ref = chars1,
    alt = chars2,
    is_mismatch = chars1 != chars2 & chars1 != "-" & chars2 != "-",
    is_gap_ref = chars1 == "-",
    is_gap_alt = chars2 == "-",
    counted = chars1 != chars2,
    stringsAsFactors = FALSE
  )
}

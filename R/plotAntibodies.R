#' Plot Antibody Data with Optional Antigen-Level Table or Time-Series Trend
#'
#' This function generates a bar plot of SAB (Single Antigen Beads) or PRA 
#' (Panel-Reactive Antibody) results from a provided data frame or a file path. 
#' It can also plot MFI values over time.
#'
#' @param result_file A data frame, a list of data frames (for trend plot), or 
#' a character string specifying the path to a file in CSV, XLS, or XLSX format.
#' @param type Character. The type of assay, either "SAB" or "PRA". Defaults to 
#' "SAB".
#' @param class Character. For PRA plots, the class of the assay, either "I" or 
#' "II". Defaults to "I".
#' @param plot_trend Logical. If TRUE, a time-series plot is generated. Defaults 
#' to FALSE.
#' @param bead_cutoffs Numeric vector. Cutoff values for categorizing MFI values.
#' Defaults to `c(2000, 1000, 500, 250)` for SAB and `c(1500, 1000, 500, 250)` 
#' for PRA.
#' @param highlight_threshold Numeric. MFI threshold for highlighting alleles in 
#' the trend plot. Defaults to 2000.
#' @param vline_dates Vector of dates. Dates to draw vertical lines on the trend 
#' plot.
#' @param add_table Logical. Whether to add the antigen-level information as a
#' table to the bottom of the bar plot. Defaults to TRUE.
#' @param x_text_angle Numeric. Angle for the antigen/allele text in the table. 
#' Defaults to 90.
#' @param palette Character. A color palette name. Defaults to "spectral".
#' @param highlight_antigen Character vector. Optional antigen(s) to highlight.
#' @param ... Additional arguments passed to the ggplot theme.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' # Plot SAB Class I data
#' plotAntibodies(deepMatchR_example[[1]], type = "SAB")
#'
#' # Plot without table
#' plotAntibodies(deepMatchR_example[[1]], type = "SAB", add_table = FALSE)
#'
#' # Plot SAB Class II data
#' plotAntibodies(deepMatchR_example[[2]], type = "SAB")
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual ylab guides theme 
#' element_blank geom_tile geom_text scale_y_discrete scale_color_manual 
#' scale_size labs geom_vline geom_line
#' @importFrom patchwork plot_layout
#' @export
plotAntibodies <- function(result_file,
                           type = "SAB",
                           class = "I",
                           plot_trend = FALSE,
                           bead_cutoffs = NULL,
                           highlight_threshold = 2000,
                           vline_dates = NULL,
                           add_table = TRUE,
                           x_text_angle = 90,
                           palette = "spectral",
                           highlight_antigen = NULL,
                           ...) {

  if (plot_trend) {
    if (is.data.frame(result_file) || !is.list(result_file) || is.null(names(result_file))) {
      stop("For trend plots, 'result_file' must be a named list of data frames.")
    }

    process_fun <- if (type == "SAB") .processSAB else .processPRA

    assay_list <- lapply(names(result_file), function(nm) {
      df <- process_fun(result_file[[nm]])
      df$sample <- nm
      df
    })
    assay_long <- data.table::rbindlist(assay_list)
    assay_long[, sample_date := as.Date(sample, format = "%m/%d/%Y")]

    alleles_to_highlight <- assay_long[, .(max_val = max(NormalValue, na.rm = TRUE)), by = allele][max_val >= highlight_threshold, allele]

    assay_long[, highlight := ifelse(allele %in% alleles_to_highlight, allele, "Other")]
    assay_long[, sample := factor(sample, levels = unique(sample))]
    assay_long[, uniqueID := paste(allele, BeadID)]

    p <- ggplot(assay_long, aes(x = sample_date, y = NormalValue, group = uniqueID, color = highlight)) +
      geom_line(alpha = 0.7) +
      scale_color_manual(values = c("Other" = "gray", setNames(.colorizer(n = length(alleles_to_highlight)), alleles_to_highlight))) +
      .themeMatchR(...) +
      labs(x = "Date", y = "MFI", color = "Allele")

    if (!is.null(vline_dates)) {
      p <- p + geom_vline(xintercept = as.Date(vline_dates), lty = 2)
    }

    return(p)
  }

  if (is.null(bead_cutoffs)) {
    bead_cutoffs <- if (type == "SAB") c(2000, 1000, 500, 250) else c(1500, 1000, 500, 250)
  }

  result0 <- if (inherits(result_file, "character")) .loadData(result_file) else result_file

  result <- if (type == "SAB") .processSAB(result0) else .processPRA(result0, class = class)
  result <- data.table::as.data.table(result)

  if (type == "SAB") {
    result[, loci := regmatches(antigen, regexpr("^[^0-9]+", antigen))]
    if (all(result$loci %in% c("A", "B", "Cw"))) {
      bw.subset <- result[!is.na(bw46)]
      if (nrow(bw.subset) > 0) {
        bw.subset[, `:=`(loci = "Bw", antigen = as.numeric(sub("[A-Za-z]+", "", bw46)))]
        result <- rbind(result, bw.subset, fill = TRUE)
      }
      result[, loci := factor(loci, levels = c("A", "B", "Bw", "Cw"))]
    } else {
      DQ.subset <- result[grepl("DQA1", allele)]
      if(nrow(DQ.subset) > 0) DQ.subset[, `:=`(loci = "DQA1", antigen = sub(".*\\*([0-9]{2}:[0-9]{2}).*", "\\1", allele))]

      DP.subset <- result[grepl("DPA1", allele)]
      if(nrow(DP.subset) > 0) DP.subset[, `:=`(loci = "DPA1", antigen = sub(".*\\*([0-9]{2}:[0-9]{2}).*", "\\1", allele))]

      result <- rbind(result, DQ.subset, DP.subset, fill = TRUE)
      result[, loci := factor(loci, levels = c("DR", "DQA1", "DQ", "DPA1", "DP"))]
    }
  } else { # PRA
    if (class == "I") {
      result[, loci := loci_family]
      result[grepl("bw", loci, ignore.case = TRUE), loci := "Bw"]
      result[, antigen := sub("Bw", "", antigen)]
      data.table::setorder(result, BeadID, loci)
      result[, group := interaction(loci, pairs)]
      
      custom_order <- c("A.1", "A.2", "B.1", "B.2", "Bw.1", "Bw.2", "C.1", "C.2")
    } else {
      result[, loci := sub("[0-9].*", "", antigen)]
      result[grepl("DR51|DR52|DR53", antigen), loci := "DR5"]

      DQ.subset <- result[grepl("DQA1", allele)]
      if(nrow(DQ.subset) > 0) DQ.subset[, `:=`(loci = "DQA1", antigen = sub(".*\\*([0-9]{2}:[0-9]{2}).*", "\\1", allele))]

      DP.subset <- result[grepl("DPA1", allele)]
      if(nrow(DP.subset) > 0) DP.subset[, `:=`(loci = "DPA1", antigen = sub(".*\\*([0-9]{2}:[0-9]{2}).*", "\\1", allele))]

      result <- rbind(result, DQ.subset, DP.subset, fill = TRUE)
      result[, group := interaction(loci, pairs)]
      custom_order <- c("DR.1", "DR.2", "DR5.1", "DR5.2", "DQ.1", "DQ.2", "DQA1.1", "DQA1.2", "DP.1", "DP.2","DPA1.1", "DPA1.2")
    }
  }

  categories <- c("Below Threshold", paste0("Level ", seq_along(bead_cutoffs)))
  result[, category := categories[findInterval(NormalValue, vec = c(-Inf, sort(bead_cutoffs), Inf), rightmost.closed = TRUE)]]

  color.palette <- .colorizer(palette = palette, n = length(bead_cutoffs) + 1)
  color.palette <- color.palette[which(rev(categories) %in% unique(result$category))]

  bead.result <- unique(result[, c("BeadID", "NormalValue", "category")])
  main_plot <- ggplot(bead.result, aes(x = reorder(BeadID, -NormalValue), y = NormalValue)) +
    geom_bar(stat = "identity", aes(fill = category), color = "black", lwd = 0.2, width = 0.7) +
    scale_fill_manual(values = rev(color.palette)) +
    .themeMatchR(...) +
    ylab("MFI Values") +
    guides(fill = "none") +
    theme(plot.background = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

  if (add_table) {
    if (!is.null(highlight_antigen)) {
      if (all(!any(grepl(paste0(result$antigen, collapse = "|"), highlight_antigen)) & !highlight_antigen %in% c("Bw4", "Bw6"))) {
        stop("highlight_antigen selection is not within the resulting data.frame")
      }
      if (all(highlight_antigen %in% c("Bw4", "Bw6"))) {
        result[, highlight := bw46 %in% highlight_antigen]
      } else {
        if (all(grepl("[*]", highlight_antigen))) {
          if (type == "PRA" && class == "II") {
            result[, highlight := allele %in% highlight_antigen & grepl(":", antigen)]
          } else {
            result[, highlight := allele %in% highlight_antigen]
          }
        } else {
          result[, highlight := antigen %in% highlight_antigen]
        }
      }
    } else {
      result[, highlight := FALSE]
    }

    if (type == "SAB") {
      table_plot <- ggplot(result, aes(x = reorder(BeadID, -NormalValue), y = loci)) +
        geom_tile(fill = "white") +
        geom_text(aes(label = antigen, color = highlight), angle = x_text_angle, size = 1.5) +
        scale_y_discrete(limits = rev) +
        .themeMatchR(..., grid_lines = "No") +
        theme(plot.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
        guides(color = "none")
    } else { # PRA
      if (class == "II") {
        result[, antigen := sub("^(DR|DQ|DP)", "", antigen)]
        ratio <- 1
      } else {
        result[, antigen := sub("^(A|B|Bw|Cw)", "", antigen)]
        ratio <- 2
      }
      result[, group := factor(group, levels = rev(custom_order))]
      result[, sizing := 2/ifelse(grepl(":", antigen), nchar(antigen), 1)]
      axis_labels <- gsub("\\..*", "", levels(result$group))

      table_plot <- ggplot(result, aes(x = reorder(BeadID, -NormalValue), y = group)) +
        geom_tile(fill = "white") +
        geom_text(aes(label = antigen, color = highlight, size = sizing), angle = x_text_angle) +
        scale_y_discrete(labels = axis_labels) +
        .themeMatchR(..., grid_lines = "No") +
        theme(plot.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
        scale_size(range = c(1.2, 2)) +
        guides(color = "none", size = "none")
    }

    combined_plot <- if (type == "SAB") main_plot / table_plot else main_plot / table_plot + plot_layout(heights = c(ratio, 1), ncol = 1)
    return(combined_plot)
  }

  return(main_plot)
}

#' @rdname plotAntibodies
#' @aliases plotSAB
#' @export
plotSAB <- function(..., type = "SAB") {
  plotAntibodies(..., type = "SAB")
}

#' @rdname plotAntibodies
#' @aliases plotPRA
#' @export
plotPRA <- function(..., type = "PRA") {
  plotAntibodies(..., type = "PRA")
}

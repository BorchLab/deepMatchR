#' Plot PRA Data with Optional Antigen-Level Table
#'
#' This function generates a bar plot of PRA (Panel-Reactive Antibody) results 
#' using a provided data frame or a file path. It processes the input data by 
#' cleaning and categorizing MFI values based on specified cutoffs and applies 
#' a chosen color palette. Optionally, it can add an antigen-level table below 
#' the plot with specific antigens highlighted.
#'
#' @param result_file A data frame containing PRA results or a character string specifying
#'   the path to a PRA file in CSV, XLS, or XLSX format.
#' @param bead_cutoffs Numeric vector. Cutoff values for categorizing MFI 
#' values. Defaults to \code{c(1500, 1000, 500, 250)}.
#' @param add_table Logical. Whether to add the antigen-level information as a 
#' table to the bottom of the plot. Defaults to \code{TRUE}.
#' @param palette Character. A color palette name (from \link[grDevices]{hcl.pals}) 
#' or a custom palette function to use for the plot. Defaults to \code{"spectral"}.
#' @param highlight_antigen Character vector. Optional antigen(s) to highlight 
#' in the table. If provided, matching antigens will be highlighted in red. 
#' Defaults to \code{NULL}.
#'
#' @return A \code{ggplot} object representing the PRA plot (and table, 
#' if \code{add_table} is \code{TRUE}).
#'
#' @examples
#' # Example using a data frame
#' plotPRA(deepMatchR_example[[1], 
#'         bead_cutoffs = c(1500, 1000, 500, 250), 
#'         add_table = TRUE, 
#'         palette = "spectral")
#'
#' @export
plotPRA <- function(result_file,
                    bead_cutoffs = c(1500, 1000, 500, 250), 
                    class = "I",
                    add_table = TRUE,
                    palette = "spectral",
                    highlight_antigen = NULL) {
  
  # 1. Read in data (data frame or file path)
  if (inherits(result_file, "character")) {
    result0 <- .loadData(result_file)
  } else {
    result0 <- result_file
  }
  
  # 2. Check and clean PRA results
  .checkSAB(result0)
  result <- .processPRA(result0, class = class)
  
  # Process loci and antigen information based on SAB type
  if (class == "I") {
    result$loci[grep("Bw", result$antigen)] <- "Bw"
    result$antigen <- str_remove(result$antigen, "Bw")
    result$group <- interaction(result$loci, result$pairs)
    custom_order <- c("A.1", "A.2", "B.1", "B.2", "Bw.1", "Bw.2", "C.1", "C.2")
  } else {
    result <- result %>%
      dplyr::mutate(loci = stringr::str_extract(antigen, "^[^0-9]+"))
    
    result$loci[grep("DR51|DR52|DR53", result$antigen)] <- "DR5"
    
    DQ.subset <- result %>%
      dplyr::filter(grepl("DQA1", allele)) %>%
      dplyr::mutate(loci = "DQA1",
                    antigen = stringr::str_extract(allele, "(?<=\\*)[0-9]{2}:[0-9]{2}"))
    
    DP.subset <- result %>%
      dplyr::filter(grepl("DPA1", allele)) %>%
      dplyr::mutate(loci = "DPA1",
                    antigen = stringr::str_extract(allele, "(?<=\\*)[0-9]{2}:[0-9]{2}"))         
    
    result <- rbind.data.frame(result, DQ.subset)
    result <- rbind.data.frame(result, DP.subset)
    result$group <- interaction(result$loci, result$pairs)
    custom_order <- c("DR.1", "DR.2", "DR5.1", "DR5.2", "DQ.1", "DQ.2", "DQA1.1", "DQA1.2", "DP.1", "DP.2", "DPA1.1", "DPA1.2")
  }
  
  # Generate categories dynamically based on bead_cutoffs, with "Below Threshold" as the first category.
  categories <- c("Below Threshold", purrr::map2_chr(bead_cutoffs, seq_along(bead_cutoffs), ~ paste0("Level ", .y)))
  
  # Categorize MFI values using findInterval to cover open-ended ranges
  result <- result %>%
    dplyr::mutate(category = categories[findInterval(NormalValue, vec = c(-Inf, sort(bead_cutoffs), Inf), rightmost.closed = TRUE)])
  
  # Generate the color palette using an internal helper function
  color.palette <- .colorizer(palette = palette, n = length(bead_cutoffs) + 1)
  
  # Prepare data for the main bar plot
  bead.result <- unique(result[, c("BeadID", "NormalValue", "category")])
  main_plot <- ggplot2::ggplot(bead.result, ggplot2::aes(x = reorder(BeadID, -NormalValue), y = NormalValue)) + 
    ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = category), 
                      color = "black", 
                      lwd = 0.2, 
                      width = 0.7) + 
    ggplot2::scale_fill_manual(values = rev(color.palette)) + 
    .dmrTheme() + 
    ggplot2::ylab("MFI Values") + 
    ggplot2::guides(fill = "none") + 
    ggplot2::theme(plot.background = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(), 
                   axis.ticks.x = ggplot2::element_blank())
  
  # Optionally add antigen-level table as an additional plot component
  if (add_table) {
    if (!is.null(highlight_antigen)) {
      if (all(!any(grepl(paste0(result$antigen, collapse = "|"), highlight_antigen)) & !highlight_antigen %in% c("Bw4", "Bw6"))) {
        stop("highlight_antigen selection is not within the resulting data.frame")
      }
      if(all(highlight_antigen %in% c("Bw4", "Bw6"))) {
        result$highlight <- result$bw46 %in% highlight_antigen
      } else {
        if(all(grepl("[*]", highlight_antigen))) {
          if(class == "II") {
            result$highlight <- result$allele %in% highlight_antigen & grepl(":", result$antigen)
          } else {
            result$highlight <- result$allele %in% highlight_antigen
          }
        } else {
          result$highlight <- result$antigen %in% highlight_antigen
        }
      }
    } else {
      result$highlight <- FALSE
    }
    
    # Removing loci prefix on antigen
    if (class == "II") {
      result$antigen <- sub("^(DR|DQ|DP)", "", result$antigen)
      ratio <- 1
    } else {
      result$antigen <- sub("^(A|B|Bw|Cw)", "", result$antigen)
      ratio <- 2
    }
    # Set order of axis labels
    result$group <- factor(result$group, levels = rev(custom_order))
    result$sizing <- 1/ifelse(grepl(":", result$antigen), nchar(result$antigen), 1)
    # Strip out the numeric portion of labels
    axis_labels <- gsub("\\..*", "", levels(result$group))
    
    table_plot <- ggplot2::ggplot(result, ggplot2::aes(x = reorder(BeadID, -NormalValue), y = group)) + 
      ggplot2::geom_tile(fill = "white") + 
      ggplot2::geom_text(ggplot2::aes(label = antigen, color = highlight, size = sizing), 
                         angle = 90) + 
      ggplot2::scale_y_discrete(labels = axis_labels) +
      .dmrTheme() + 
      ggplot2::theme(plot.background = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(), 
                     axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + 
      scale_size(range = c(1.2, 2)) + 
      ggplot2::guides(color = "none", size = "none")
    
    # Combine the main plot and the table plot using patchwork 
    combined_plot <- main_plot / table_plot + plot_layout(heights = c(ratio, 1), ncol = 1)
    return(combined_plot)
  }
  
  return(main_plot)
}

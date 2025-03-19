#' @param result_file A data frame containing SAB results or a character string specifying
#'   the path to a SAB file in CSV, XLS, or XLSX format.
#' @param bead_cutoffs Numeric vector. Categories of colors to use for visualizing 
#' groups based on the MFI values.
#' @param add_table Logical. Whether to add the antigen-level information in the form
#' of a table to the bottom on the plot. 
#' @param palette Character. A color palette name (from \link[grDevices]{hcl.pals}) or a custom
#'   palette function to use for the plot. Defaults to \code{"spectral"}.

plotSAB <- function(result_file,
                    bead_cutoffs = c(2000, 1000, 500, 250), 
                    add_table = TRUE
                    palette = "spectral")
# 1. Read in data (data frame or file path)
  if (inherits(result_file, "character")) {
    result0 <- .loadData(result_file)
  } else {
    result0 <- result_file
  }
  
  # 2. Check and Clean SAB results
  .checkSAB(result0)
  result <- .processSAB(result0)
  
  if(all(grep("A|B|C", result$loci))) {
    result <- result %>%
      mutate(loci = str_extract(allele, "^[^*]+"))
    bw.subset <- result %>%
                  subset(!is.na(bw46)) %>%
                  mutate(loci = "Bw",
                         antigen = as.numeric(sub("[A-Za-z]+", "", bw46)))
    result <- rbind.data.frame(result, bw.subset)
    result$loci <- factor(result$loci, levels = c("A", "B", "Bw", "C"))
  } else {
    result <- result %>%
      mutate(loci = str_extract(antigen, "^[^0-9]+"))
    
    DQ.subset <- result %>%
      subset(grepl("DQA1", allele)) %>%
      mutate(loci = "DQA1",
             antigen = str_extract(allele, "(?<=\\*)[0-9]{2}:[0-9]{2}"))
    
    DP.subset <- result %>%
      subset(grepl("DPA1", allele)) %>%
      mutate(loci = "DPA1",
             antigen = str_extract(allele, "(?<=\\*)[0-9]{2}:[0-9]{2}"))         
            
    result <- rbind.data.frame(result, DQ.subset)
    result <- rbind.data.frame(result, DP.subset)
    result$loci <- factor(result$loci, levels = c("DR", "DQA1", "DQ", "DPA1", "DP"))
  }
  
  
  
  #Generate categories dynamically based on cut.offs, with Above Threshold at the end
  categories <- c("Below Threshold", map2_chr(bead_cutoffs, seq_along(bead_cutoffs), ~ paste0("Level ", .y)))
  
  # Adjust cut.offs to include -Inf and Inf to cover the open-ended ranges
  result <- result %>%
    mutate(
      category = categories[findInterval(NormalValue, vec = c(-Inf, sort(bead_cutoffs), Inf), rightmost.closed = TRUE)]
    )
  
  # Generate the color palette using internal helper function
  color.palette <- .colorizer(palette = palette, 
                              n = length(unique(result$category)))
  
  bead.result <- unique(result[,c("BeadID", "NormalValue", "category")])
  plot <- ggplot(bead.result, aes(x = reorder(BeadID, -NormalValue), y = NormalValue)) + 
    geom_bar(stat = "identity", aes(fill = category), 
             color = "black", 
             lwd = 0.2, 
             width = 0.7) + 
    scale_fill_manual(values = rev(color.palette)) + 
    theme_clean() + 
    ylab("MFI Values") + 
    guides(fill = "none") + 
    theme(plot.background = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
  
  if(add_table) {
    if (!is.null(highlight_antigen)) {
      if (any(grepl(highlight_antigen, result$antigen))) {
        result$highlight <- apply(result, 1, function(row) {
          any(highlight_antigen == row["antigen"])
        })
      } else {
        stop("highlight_antigen selection is not within the data.frame")
      }
    } else {
      table$highlight <- FALSE
    }
    
    table.plot <- ggplot(result, aes(x = reorder(BeadID, -NormalValue), y = loci)) + 
      geom_tile(fill = "white") + 
      geom_text(aes(label = antigen, color = highlight), 
                angle = 90, 
                size = 1.5) + 
      scale_y_discrete(limits=rev) + 
      theme_clean() + 
      theme(plot.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank()) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) + 
      guides(color = "none")
    plot <- plot + table.plot
  }
  return(plot)
}

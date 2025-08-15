# Basic theme for all plots
.themeMatchR <- function(base_size = 12,
                         base_family = "sans",
                         grid_lines = "Y",
                         axis_lines = FALSE,
                         legend_position = "right") {
  
  t <- ggplot2::theme_bw(base_size = base_size, base_family = base_family)
  t <- t %+replace%
    ggplot2::theme(
      # Plot titles and caption
      plot.title = ggplot2::element_text(
        size = rel(1.2), hjust = 0, face = "bold",
        margin = ggplot2::margin(b = base_size / 2)
      ),
      plot.subtitle = ggplot2::element_text(
        size = rel(1.0), hjust = 0,
        margin = ggplot2::margin(b = base_size)
      ),
      plot.caption = ggplot2::element_text(
        size = rel(0.8), hjust = 1, color = "grey50",
        margin = ggplot2::margin(t = base_size / 2)
      ),
      
      # Backgrounds and borders
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.75),
      
      # Remove all grid lines by default; they will be added back conditionally
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      
      # Axis text, titles, and ticks
      axis.title = ggplot2::element_text(size = rel(1.0)),
      axis.text = ggplot2::element_text(size = rel(0.9), color = "black"),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.5),
      
      # Legend customization
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = rel(0.9), face = "bold"),
      legend.text = ggplot2::element_text(size = rel(0.85)),
      legend.position = legend_position,
      
      # Facet (strip) customization
      strip.background = ggplot2::element_rect(fill = "grey90", color = "black", linewidth = 0.75),
      strip.text = ggplot2::element_text(
        size = rel(1.0), face = "bold", color = "black",
        margin = ggplot2::margin(t = base_size / 4, b = base_size / 4)
      )
    )
  
  # Conditionally add major grid lines based on the 'grid_lines' parameter
  grid_lines <- toupper(grid_lines)
  if (grid_lines %in% c("Y", "XY")) {
    t <- t + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "grey85", linewidth = 0.5))
  }
  if (grid_lines %in% c("X", "XY")) {
    t <- t + ggplot2::theme(panel.grid.major.x = ggplot2::element_line(color = "grey85", linewidth = 0.5))
  }
  
  # Conditionally add axis lines
  if (axis_lines) {
    t <- t + ggplot2::theme(axis.line = ggplot2::element_line(color = "black", linewidth = 0.5))
  }
  
  return(t)
}

# Check if nessecary columns are present in SAB results
.checkSAB <- function(file) {
  if(!all(c("BeadID", "SpecAbbr", "Specificity", "NormalValue") %in% colnames(file))) {
    stop("Please ensure the SAB file includes the 'BeadID', 'SpecAbbr', 'Specificity', 'NormalValue'")
  }
}

# Loading csv or xls/xlsx data
#' @importFrom readxl read_excel
.loadData <- function(file_path) {
  # Ensure the file exists
  if (!file.exists(file_path)) {
    stop("The specified file does not exist.")
  }
  
  # Extract file extension
  file_ext <- tools::file_ext(file_path)
  
  # Load the appropriate library
  if (file_ext == "csv") {
    return(read.csv(file_path, stringsAsFactors = FALSE))
  } else if (file_ext %in% c("xls", "xlsx")) {
    return(readxl::read_excel(file_path))
  } else {
    stop("Unsupported file format. Please provide a .csv, .xls, or .xlsx file.")
  }
}

.strExtract <- function(string, pattern) {
  # regexpr finds the match position and length
  match_info <- regexpr(pattern, string)
  
  # Extract the substring using the start position and length
  # substr preserves vector length, returning "" for non-matches (where match_info is -1)
  extracted <- substr(
    string,
    match_info,
    match_info + attr(match_info, "match.length") - 1
  )
  
  # Replace the empty strings from non-matches with NA
  ifelse(extracted == "", NA_character_, extracted)
}

#Pulling a color palette for visualizations
#' @importFrom grDevices hcl.colors
#' @keywords internal
.colorizer <- function(palette = "spectral", 
                       n= NULL) {
  colors <- hcl.colors(n=n, palette = palette, fixup = TRUE)
  return(colors)
}

#' @importFrom dplyr filter mutate select distinct arrange group_by ungroup 
#'   summarise relocate left_join n
#' @importFrom tidyr unnest
.processSAB <- function(result0) {
  result <- result0 %>%
    dplyr::select(BeadID, SpecAbbr, Specificity, NormalValue) %>%
    dplyr::distinct(Specificity, .keep_all = TRUE) %>%
    mutate(
      antigen = .strExtract(SpecAbbr, '[ABCDRQP][[:alnum:]]+'),
      bw46 = .strExtract(SpecAbbr, 'Bw[46]'),
      Specificity_truncated = .strExtract(Specificity, '[ABCD][^ ()]*[0-9]')
    ) %>%
    mutate(
      allele = strsplit(gsub(",-,", "_", Specificity_truncated, fixed = TRUE), "_")
    ) %>%
    dplyr::select(-SpecAbbr, -Specificity, -Specificity_truncated) %>%
    relocate(BeadID, antigen, bw46, allele, NormalValue) %>%
    tidyr::unnest(allele) %>%
    mutate(loci = sub("\\*.*", "", allele)) %>%
    mutate(mfi_min = min(NormalValue), .by = allele) %>%
    arrange(allele, desc(NormalValue)) %>%
    filter(!is.na(allele))
  return(result)
}

# Helper function to replicate stringr::str_pad(side = "right")
.padRightBase <- function(vec, len, pad = "-") {
  pad_lengths <- pmax(0, len - nchar(vec))
  padding <- sapply(pad_lengths, function(n) paste(rep(pad, n), collapse = ""))
  paste0(vec, padding)
}

#' @importFrom data.table as.data.table `:=`
.processPRA <- function(result0, class = "I") {
  #Remove Empty Cells
  result0 <- result0[!result0$SpecAbbr == "",]
  
  # Identify class I vs II using first token of SpecAbbr
  Spec.pattern <- do.call(rbind, strsplit(result0$SpecAbbr, ",", fixed = TRUE))
  classI.pos <- grep("A", Spec.pattern[, 1])
  
  if (class == "I") {
    result <- result0[classI.pos, ]
  } else {
    result <- result0[-classI.pos, ]
  }
  
  dt <- as.data.table(result)[, .(BeadID, SpecAbbr, Specificity, NormalValue)]
  dt[, antigen_vec := strsplit(SpecAbbr, ",", fixed = TRUE)]
  dt[, allele_vec  := strsplit(Specificity, ",", fixed = TRUE)]
  dt[, rid := .I]  # stable row id per bead row
  
  # Expand one row at a time, padding shorter side with "-"
  expanded <- dt[, {
    ant <- antigen_vec[[1]]
    all <- allele_vec[[1]]
    max_len <- max(length(ant), length(all))
    if (length(ant) < max_len) ant <- c(ant, rep("-", max_len - length(ant)))
    if (length(all) < max_len) all <- c(all, rep("-", max_len - length(all)))
    data.table(position = seq_len(max_len), antigen = ant, allele = all)
  }, by = .(BeadID, NormalValue, rid)]
  
  # Class-specific position logic (mirrors your dplyr version)
  if (class != "I") {
    # Duplicate rows with position >= 5, shifting by +2
    dup <- expanded[position >= 5][, position := position + 2]
    # Shift positions 7 and 8 by +4 in the original set
    expanded[position %in% c(7, 8), position := position + 4]
    expanded <- rbindlist(list(expanded, dup), use.names = TRUE)
    string.pattern <- "-"   # placeholders are "-" in this representation
  } else {
    # For Class I, positions > 4 are shifted down by 2
    expanded[position > 4, position := position - 2]
    string.pattern <- "-"
  }
  
  # Impute antigen when placeholder but allele is not, using previous antigen within bead row
  setorder(expanded, BeadID, rid, position)
  expanded[, prev_antigen := shift(antigen), by = .(BeadID, rid)]
  expanded[antigen == string.pattern & allele != string.pattern,
           antigen := prev_antigen]
  expanded[, prev_antigen := NULL]
  
  # Clean and annotate
  out <- expanded[
    antigen != string.pattern & allele != string.pattern & antigen != "" & allele != "",
    .(BeadID, antigen, allele, NormalValue)
  ]
  
  out[, antigen := gsub("-", "", antigen)]
  out[, allele  := gsub("-", "", allele)]
  out[, bw46    := ifelse(grepl("Bw[46]", antigen), antigen, NA_character_)]
  out[, loci    := sub("\\*.*", "", allele)]
  out[, mfi_min := min(NormalValue, na.rm = TRUE), by = allele]
  
  # Pair flag and de-dup
  setorder(out, BeadID, loci)
  out[, pairs := rep(c(1, 2), length.out = .N), by = BeadID]
  out <- unique(out, by = c("BeadID", "antigen", "allele"))
  
  # Final column order
  setcolorder(out, c("BeadID", "antigen", "bw46", "allele", "loci", "NormalValue", "pairs"))
  return(as.data.frame(out))
}


.alphanumericalSort <- function(x, ignore.case = TRUE) {
  if (!is.character(x)) {
    message("Input 'x' is not a character vector. Attempting to convert.")
    x <- as.character(x)
  }
  unique_elements <- unique(x)
  alpha_part <- gsub("[0-9]+", "", unique_elements, perl = TRUE)
  if (ignore.case) {
    alpha_part <- tolower(alpha_part)
  }
  numeric_part <- suppressWarnings(as.numeric(gsub("[^0-9]", "", unique_elements)))
  sorted_levels <- unique_elements[
    order(alpha_part, numeric_part)
  ]
  
  return(sorted_levels)
}

#' @importFrom keras3 load_model
.loadModel <- function(chain, class) {
  select  <- system.file("extdata", paste0(class, "_encoder.keras"), 
                         package = "deepMatchR")
  model <- load_model(select, compile = FALSE)
  return(model)
}

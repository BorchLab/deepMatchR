#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @importFrom dplyr summarize row_number
#' @importFrom ggplot2 %+replace% rel coord_flip
#' @importFrom data.table data.table set setorder fifelse rbindlist setcolorder
#' @importFrom treemapify geom_treemap geom_treemap_text geom_treemap_subgroup_border geom_treemap_subgroup_text
#' @importFrom utils head read.csv globalVariables
NULL

# Quiet R CMD check notes about non-standard evaluation
if(getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(".", "BeadID", "SpecAbbr", "Specificity", "NormalValue", "antigen_vec",
      "allele_vec", "rid", ".I", "position", "is_bw", "antigen", "is_cant",
      "allele", "prev_antigen", "shift", "allele_locus", "loci_family0",
      "loci_family", "bw_label", "bw46", "mfi_min", "pairs", ".N",
      "Specificity_truncated", "desc", "count", "positive_count", "subtotal",
      "percent_positive", "pp_max", "loci", "AUC", "norm_AUC", ".data",
      "deepMatchR_cregs", "creg", "deepMatchR_eplets", "eplet", "median",
      "max_val", "sample_date", "highlight", "setNames", "reorder", "category",
      "group", "sizing", "positive.bead", "count_above", "count_total",
      "row_number", "sym")
  )
}

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

#TODO fix DPA reporting
#' @importFrom data.table as.data.table `:=`
.processPRA <- function(result0, class = "I") {
  
  #Remove Empty Cells
  result0 <- result0[!result0$SpecAbbr == "",]
  
  # --- Split incoming table into Class I vs II based on first SpecAbbr token ---
  Spec.pattern <- do.call(rbind, strsplit(result0$SpecAbbr, ",", fixed = TRUE))
  classI.pos   <- grep("A", Spec.pattern[, 1])  # "A" for HLA-A at first slot
  
  if (class == "I") {
    result <- result0[classI.pos, ]
    string.pattern <- "--------"    # original placeholder for Class I
  } else {
    result <- result0[-classI.pos, ]
    string.pattern <- "------------" # original placeholder for Class II
  }
  
  # --- Expand rows into positions, padding with class-specific placeholder ---
  dt <- as.data.table(result)[, .(BeadID, SpecAbbr, Specificity, NormalValue)]
  dt[, antigen_vec := strsplit(SpecAbbr, ",", fixed = TRUE)]
  dt[, allele_vec  := strsplit(Specificity, ",", fixed = TRUE)]
  dt[, rid := .I]  # stable row id per bead
  
  expanded <- dt[, {
    ant <- trimws(antigen_vec[[1]])
    all <- trimws(allele_vec[[1]])
    max_len <- max(length(ant), length(all))
    if (length(ant) < max_len) ant <- c(ant, rep(string.pattern, max_len - length(ant)))
    if (length(all) < max_len) all <- c(all, rep(string.pattern, max_len - length(all)))
    
    # If antigen is a literal "-", mark BOTH sides as placeholder so the row drops later
    idx_dash <- ant == "-"
    if (any(idx_dash)) {
      ant[idx_dash] <- string.pattern
      all[idx_dash] <- string.pattern
    }
    
    data.table(position = seq_len(max_len),
               antigen  = ant,
               allele   = all)
  }, by = .(BeadID, NormalValue, rid)]
  
  # --- Class-specific position rules (mirror original dplyr logic) ---
  if (class != "I") {
    dup <- expanded[position >= 5][, position := position + 2]
    expanded[position %in% c(7, 8), position := position + 4]
    expanded <- rbindlist(list(expanded, dup), use.names = TRUE)
  } else {
    expanded[position > 4, position := position - 2]
  }
  
  # --- Class I only: move C* alleles from Bw rows onto C antigen rows (by order), then clear Bw alleles
  if (class == "I") {
    expanded[, is_bw := grepl("^Bw[46]$", antigen)]
    expanded[, is_cant := grepl("^(Cw|C)", antigen)]
    
    expanded[, {
      # indexes local to this group
      rows <- .I
      bw_idx <- which(is_bw & grepl("^C\\*", allele))
      c_idx  <- which(is_cant & allele == string.pattern)
      if (length(bw_idx) && length(c_idx)) {
        bw_idx <- bw_idx[order(position[bw_idx])]
        c_idx  <- c_idx[order(position[c_idx])]
        k <- min(length(bw_idx), length(c_idx))
        # move alleles from Bw -> C; clear Bw alleles
        set(expanded, i = rows[c_idx[seq_len(k)]], j = "allele", value = expanded$allele[rows[bw_idx[seq_len(k)]]])
        set(expanded, i = rows[bw_idx[seq_len(k)]], j = "allele", value = string.pattern)
      }
      NULL
    }, by = .(BeadID, rid)]
  }
  
  # --- Impute antigen ONLY when placeholder is present (not "-") ---
  setorder(expanded, BeadID, rid, position)
  expanded[, prev_antigen := shift(antigen), by = .(BeadID, rid)]
  expanded[antigen == string.pattern & allele != string.pattern, antigen := prev_antigen]
  expanded[, prev_antigen := NULL]
  
  # --- Keep valid rows; strip hyphens and trim; keep rid to help synthesize BW later ---
  out <- expanded[
    antigen != string.pattern & allele != string.pattern &
      antigen != "-" & allele != "-" &
      antigen != ""  & allele != "",
    .(BeadID,
      rid,
      antigen = trimws(gsub("-", "", antigen)),
      allele  = trimws(gsub("-", "", allele)),
      NormalValue)
  ]
  
  # --- Locus annotations ---
  out[, allele_locus := sub("\\*.*", "", allele)]                # e.g., DRB1, DQA1, B, C, etc.
  out[, loci_family0 := toupper(sub("[0-9].*", "", antigen))]    # e.g., DR, DQ, DP, A, B, Cw, Bw...
  # normalize "CW" -> "C", keep "BW" as "BW"
  out[, loci_family := fifelse(loci_family0 == "CW", "C", loci_family0)]
  out[, loci_family0 := NULL]
  
  # --- Class II guard-rails + DR5 mapping ---
  if (class != "I") {
    # Map DR52/DR53 (DRB3/4/5) into DR5 family
    is_dr5_antigen <- grepl("^DR(51|52|53)$", toupper(out[,antigen])) | grepl("^DR5(1|2|3)$", toupper(out[,antigen]))
    is_dr5_allele  <- grepl("^DRB[345]$", out[,allele_locus])
    out[is_dr5_antigen | is_dr5_allele, loci_family := "DR5"]
    
    # Guard-rails
    out[loci_family == "DQ"  & !grepl("^(DQA1|DQB1)$", allele_locus),
        c("allele","allele_locus") := .(NA_character_, NA_character_)]
    out[loci_family == "DP"  & !grepl("^(DPA1|DPB1)$", allele_locus),
        c("allele","allele_locus") := .(NA_character_, NA_character_)]
    out[loci_family %in% c("DR","DR5") & !grepl("^DRB", allele_locus),
        c("allele","allele_locus") := .(NA_character_, NA_character_)]
    out <- out[!is.na(allele) & allele != ""]
  }
  
  # --- BW synthesis for Class I ---
  if (class == "I") {
    # Collect Bw labels present per bead from expanded (post-shift)
    bw_map <- unique(expanded[grepl("^Bw[46]$", antigen), .(BeadID, rid, bw_label = antigen)])
    # Remove vendor Bw rows (they had their C* alleles already moved off)
    out <- out[!grepl("^Bw[46]$", antigen)]
    # For each bead/rid with a Bw label, add a clean BW row (do not touch A/B/C)
    if (nrow(bw_map)) {
      # Build one BW row per label using the bead's NormalValue (take max per bead/rid for safety)
      nv_map <- out[, .(NormalValue = max(NormalValue, na.rm = TRUE)), by = .(BeadID, rid)]
      bw_synth <- merge(bw_map, nv_map, by = c("BeadID","rid"), all.x = TRUE)
      bw_synth[, `:=`(
        antigen      = bw_label,       # "Bw4" / "Bw6"
        allele       = bw_label,       # simple label for plotting
        allele_locus = "BW",
        loci_family  = "BW",
        bw46         = bw_label
      )]
      bw_synth <- bw_synth[, .(BeadID, rid, antigen, allele, NormalValue, allele_locus, loci_family, bw46)]
      # rbind to out (align columns)
      missing_cols <- setdiff(names(out), names(bw_synth))
      if (length(missing_cols)) bw_synth[, (missing_cols) := NA]
      out <- rbindlist(list(out, bw_synth), use.names = TRUE, fill = TRUE)
    }
  }
  
  # --- Metrics ---
  out[, bw46 := fifelse(grepl("^Bw[46]$", antigen), antigen, NA_character_)]
  out[, mfi_min := min(NormalValue, na.rm = TRUE), by = allele]
  
  # --- Pairs assignment ---
  if (class == "I") {
    setorder(out, BeadID, loci_family, antigen, allele)
    out[, pairs := rep(c(1, 2), length.out = .N), by = .(BeadID, loci_family)]
  } else {
    setorder(out, BeadID, loci_family, antigen, allele_locus)
    out[, pairs := as.integer(factor(antigen, levels = unique(antigen))), by = .(BeadID, loci_family)]
  }
  
  # De-dup like dplyr::distinct(BeadID, antigen, allele, .keep_all = TRUE)
  out <- unique(out, by = c("BeadID", "antigen", "allele"))
  
  # Final column order; drop helper 'rid'
  setcolorder(out, c("BeadID","antigen","bw46","allele","allele_locus","loci_family","NormalValue","mfi_min","pairs","rid"))
  out[, rid := NULL]
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

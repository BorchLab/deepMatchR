"%!in%" <- Negate("%in%")

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
#' @importFrom stringr str_extrext function str_replace_all
#' @importFrom tidyr separate_longer_delim
.processSAB <- function(result0) {
  result <- result0 %>%
    dplyr::select(BeadID, SpecAbbr, Specificity, NormalValue) %>%
    dplyr::distinct(Specificity, .keep_all = TRUE) %>%
    mutate(
      antigen               = str_extract(SpecAbbr, '[ABCDRQP][:alnum:]+'),
      bw46                  = str_extract(SpecAbbr, 'Bw[46]'),
      Specificity_truncated = str_extract(Specificity, '[ABCD].*[0-9]')
    ) %>%
    mutate(
      allele = str_replace_all(Specificity_truncated, ",-,", "_"),
      loci = str_extract(allele, "^[^*]+")
    ) %>%
    dplyr::select(-SpecAbbr, -Specificity, -Specificity_truncated) %>%
    relocate(BeadID, antigen, bw46, allele, NormalValue) %>%
    separate_longer_delim(allele, "_") %>%
    mutate(mfi_min = min(NormalValue), .by = allele) %>%
    arrange(allele, desc(NormalValue)) %>%
    filter(!is.na(allele)) 
  return(result)
}


.processPRA <- function(result0, class = "I") {
  # Step 0: Determine class-specific rows
  Spec.pattern <- str_split(result0$SpecAbbr, ",", simplify = TRUE)
  classI.pos <- grep("A", Spec.pattern[,1])
  
  if (class == "I") {
    result <- result0[classI.pos,]
  } else {
    result <- result0[-classI.pos,]
  }
  
  # Step 1: Split and pad antigen & allele vectors
  result_expanded <- result %>%
    dplyr::select(BeadID, SpecAbbr, Specificity, NormalValue) %>%
    rowwise() %>%
    mutate(
      antigen_vec = str_split(SpecAbbr, ","),
      allele_vec  = str_split(Specificity, ","),
      max_len     = max(length(antigen_vec), length(allele_vec)),
      antigen_vec = list(str_pad(antigen_vec, max_len, side = "right", pad = "-")),
      allele_vec  = list(str_pad(allele_vec,  max_len, side = "right", pad = "-"))
    ) %>%
    ungroup() %>%
    dplyr::select(BeadID, NormalValue, antigen_vec, allele_vec) %>%
    unnest_longer(antigen_vec, values_to = "antigen", indices_to = "position")
  
  # Step 2: Adjust position and duplicate if needed (Class II logic)
  if (class != "I") {
    # Duplicate rows for positions >=5 with adjusted positions
    result_expanded2 <- result_expanded %>% filter(position >= 5)
    result_expanded  <- result_expanded %>%
      mutate(position = ifelse(position %in% c(7, 8), position + 4, position))
    result_expanded2 <- result_expanded2 %>%
      mutate(position = position + 2)
    result_expanded <- bind_rows(result_expanded, result_expanded2)
    string.pattern <- "------------"
  } else {
    # For Class I, normalize position for downstream indexing
    result_expanded <- result_expanded %>%
      group_by(BeadID) %>%
      mutate(position = ifelse(position > 4, position - 2, position)) %>%
      ungroup()
    string.pattern <- "--------"
  }
  
  # Imputing antigen level info for lazy load
  result_expanded <- result_expanded %>%
    mutate(
      prev_antigen = lag(antigen)
    ) %>%
    rowwise() %>%
    mutate(
      antigen = if (antigen == string.pattern && allele_vec[[position]] != string.pattern) {
        prev_antigen
      } else {
        antigen
      }
    ) %>%
    dplyr::select(-prev_antigen) %>%
    ungroup()
  
  # Step 3: Map alleles by position and clean up
  result_expanded <- result_expanded %>%
    rowwise() %>%
    mutate(allele = allele_vec[position]) %>%
    ungroup() %>%
    mutate(pairs   = rep(c(1,2), length.out = n())) %>%
    filter(antigen != "-", allele != "-") %>%
    mutate(
      antigen = str_remove_all(antigen, "-"),
      allele  = str_remove_all(allele, "-"),
      bw46    = ifelse(str_detect(antigen, "Bw[46]"), antigen, NA_character_),
      loci    = str_extract(allele, "^[^*]+"),
      mfi_min = min(NormalValue, na.rm = TRUE), .by = allele) %>%
    filter(antigen != string.pattern, allele != "", antigen != "") %>%
    distinct(BeadID, antigen, allele, .keep_all = TRUE) %>%
    dplyr::select(BeadID, antigen, bw46, allele, loci, NormalValue, pairs)
  
  return(result_expanded)
}

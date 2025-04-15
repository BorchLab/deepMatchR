"%!in%" <- Negate("%in%")

# Basic theme for ggplot modified from plotthis R package
.dmrTheme <- function(aspect.ratio = NULL, base_size = 12, font_family = NULL, ...) {
  text_size_scale <- base_size / 12
  args1 <- list(
    aspect.ratio = aspect.ratio,
    text = element_text(size = 12 * text_size_scale, family = font_family, color = "black"),
    plot.title = element_text(size = 14 * text_size_scale, family = font_family, colour = "black", vjust = 1),
    plot.subtitle = element_text(size = 13 * text_size_scale, family = font_family, hjust = 0, margin = margin(b = 3)),
    plot.background = element_rect(fill = "white", color = "white"),
    plot.margin = margin(10, 10, 10, 10),
    axis.line = element_blank(),
    axis.title = element_text(size = 13 * text_size_scale, family = font_family, colour = "black"),
    axis.text = element_text(size = 12 * text_size_scale, family = font_family, colour = "black"),
    strip.text = element_text(size = 12.5 * text_size_scale, family = font_family, colour = "black", hjust = 0.5, margin = margin(3, 3, 3, 3)),
    strip.background = element_rect(fill = "transparent", linetype = 0),
    strip.switch.pad.grid = unit(-1, "pt"),
    strip.switch.pad.wrap = unit(-1, "pt"),
    strip.placement = "outside",
    legend.title = element_text(size = 12 * text_size_scale, family = font_family, colour = "black", hjust = 0),
    legend.text = element_text(size = 11 * text_size_scale, family = font_family, colour = "black"),
    legend.key = element_rect(fill = "transparent", color = "transparent"),
    legend.key.size = unit(10, "pt"),
    legend.background = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill = "transparent", colour = "black", linewidth = 1)
  )
  args2 <- as.list(match.call())[-1]
  call.envir <- parent.frame(1)
  args2 <- lapply(args2, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }
  args <- args1[names(args1) %in% formalArgs(theme)]
  do.call(what = theme, args = args)
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

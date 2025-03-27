"%!in%" <- Negate("%in%")

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
    filter(!is.na(allele)) %>%
    distinct(allele, .keep_all = TRUE) 
  return(result)
}

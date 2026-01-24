#' Convert HLA Alleles to Serological Equivalents
#'
#' Converts HLA alleles in IMGT/HLA nomenclature (e.g., "A*01:01") to their
#' serological equivalents (e.g., "A1") using WMDA standard nomenclature.
#'
#' @param x Character vector of HLA alleles OR an \code{hla_genotype} object.
#' @param locus Optional character. Override locus detection for ambiguous cases.
#' @param resolve_splits Logical. If TRUE (default), attempts to infer split antigens
#'   from broad antigens (e.g., DR2 -> DR15 or DR16 based on the allele).
#' @param return Character. Output format:
#'   \itemize{
#'     \item \code{"serology"} (default): Returns character vector of serology strings
#'     \item \code{"genotype"}: Returns hla_genotype with \code{*_ser_*} columns added
#'     \item \code{"data.frame"}: Returns data.frame with allele, locus, serology columns
#'   }
#' @param na_action Character. How to handle alleles without serology mapping:
#'   \itemize{
#'     \item \code{"NA"} (default): Return NA for unknown alleles
#'     \item \code{"warn"}: Return NA and issue a warning
#'     \item \code{"error"}: Stop with an error
#'   }
#'
#' @return Depends on \code{return} parameter. See above.
#'
#' @details
#' This function uses the WMDA (World Marrow Donor Association) nomenclature files
#' to map HLA alleles to their serological equivalents. The mapping considers:
#' \itemize{
#'   \item Unambiguous assignments (highest priority)
#'   \item Possible assignments
#'   \item Assumed assignments
#'   \item Expert assignments
#' }
#'
#' When \code{resolve_splits = TRUE}, broad antigens like DR2 are resolved to their
#' split antigens (DR15 or DR16) based on the specific allele.
#'
#' P-group notation (e.g., "A*01:01P") is automatically resolved to the reference allele.
#'
#' @examples
#' # Single allele conversion
#' toSerology("A*01:01")
#'
#' # Vector of alleles
#' toSerology(c("A*01:01", "B*07:02", "DRB1*03:01"))
#'
#' # Get detailed data.frame
#' toSerology(c("A*01:01", "B*07:02"), return = "data.frame")
#'
#' # With hla_genotype object
#' geno <- hlaGeno(data.frame(
#'   A_1 = "A*01:01", A_2 = "A*02:01",
#'   B_1 = "B*07:02", B_2 = "B*08:01",
#'   stringsAsFactors = FALSE
#' ))
#' toSerology(geno, return = "serology")
#' toSerology(geno, return = "genotype")
#'
#' @export
toSerology <- function(x,
                       locus = NULL,
                       resolve_splits = TRUE,
                       return = c("serology", "genotype", "data.frame"),
                       na_action = c("NA", "warn", "error")) {

  return <- match.arg(return)
  na_action <- match.arg(na_action)

  # Handle hla_genotype input
  if (inherits(x, "hla_genotype")) {
    return(.toSerologyGenotype(x, locus, resolve_splits, return, na_action))
  }

  # Validate character input
  if (!is.character(x)) {
    stop("'x' must be a character vector or hla_genotype object")
  }

  # Load WMDA data
  wmda_data <- .loadWmdaData()

  # Process each allele
  results <- lapply(x, function(allele) {
    .convertSingleAllele(allele, locus, resolve_splits, na_action, wmda_data)
  })

  results_df <- data.table::rbindlist(results)

  # Return based on format
  switch(return,
         "serology" = results_df$serology_full,
         "data.frame" = as.data.frame(results_df),
         stop("Invalid return type"))
}

#' @noRd
.toSerologyGenotype <- function(geno, locus, resolve_splits, return, na_action) {
  validateHlaGeno(geno)

  # Get all allele columns
  df <- geno$data
  allele_cols <- grep("^[A-Z]+[0-9]*_[12]$", names(df), value = TRUE)

  # Extract all alleles
  all_alleles <- unlist(df[, allele_cols, drop = FALSE], use.names = FALSE)
  all_alleles <- all_alleles[!is.na(all_alleles) & nzchar(all_alleles)]
  unique_alleles <- unique(all_alleles)

  # Convert all unique alleles
  wmda_data <- .loadWmdaData()
  conversion_list <- lapply(unique_alleles, function(allele) {
    .convertSingleAllele(allele, locus, resolve_splits, na_action, wmda_data)
  })
  conversion_df <- data.table::rbindlist(conversion_list)

  # Create lookup
  serology_lookup <- stats::setNames(conversion_df$serology_full, conversion_df$allele)

  if (return == "serology") {
    # Return just the serology values for all alleles in the genotype
    return(serology_lookup[all_alleles])
  }

  if (return == "data.frame") {
    return(as.data.frame(conversion_df))
  }

  if (return == "genotype") {
    # Add *_ser_* columns to the data
    new_df <- df
    for (col in allele_cols) {
      ser_col <- sub("_([12])$", "_ser_\\1", col)
      new_df[[ser_col]] <- serology_lookup[df[[col]]]
    }

    # Create new hla_genotype
    structure(
      list(
        data = new_df,
        locus_present = geno$locus_present
      ),
      class = "hla_genotype"
    )
  }
}

#' @noRd
.convertSingleAllele <- function(allele, locus_override, resolve_splits, na_action, wmda_data) {
  # Handle NA/empty
  if (is.na(allele) || !nzchar(trimws(allele))) {
    return(data.table::data.table(
      allele = allele,
      locus = NA_character_,
      allele_2f = NA_character_,
      serology = NA_character_,
      serology_full = NA_character_
    ))
  }

  # Parse the allele
  parsed <- .parseAllele(allele)
  if (is.null(parsed)) {
    return(.handleUnknown(allele, na_action, "Could not parse allele format"))
  }

  allele_locus <- parsed$locus
  allele_2f <- parsed$allele_2f

  # Override locus if specified
  if (!is.null(locus_override)) {
    allele_locus <- locus_override
  }

  # Resolve P-group if present
  if (grepl("P$", allele_2f)) {
    resolved <- .resolvePGroup(allele_locus, allele_2f, wmda_data$pgroups)
    if (!is.null(resolved)) {
      allele_2f <- resolved
    }
  }

  # Look up serology
  serology <- .lookupSerology(allele_locus, allele_2f, wmda_data$serology)

  if (is.na(serology)) {
    return(.handleUnknown(allele, na_action, paste0("No serology mapping for ", allele_locus, "*", allele_2f)))
  }

  # Resolve splits if requested
  if (resolve_splits) {
    serology <- .resolveBroadToSplit(allele_locus, allele_2f, serology, wmda_data)
  }

  # Get serology prefix for full name
  prefix <- .getSerologyPrefix(allele_locus)
  serology_full <- paste0(prefix, serology)

  data.table::data.table(
    allele = allele,
    locus = allele_locus,
    allele_2f = allele_2f,
    serology = serology,
    serology_full = serology_full
  )
}

#' @noRd
.parseAllele <- function(allele) {
  allele <- trimws(allele)

  # Handle format: A*01:01 or A*01:01:01 or A*01:01:01:01
  if (!grepl("\\*", allele)) {
    return(NULL)
  }

  parts <- strsplit(allele, "\\*")[[1]]
  if (length(parts) != 2) {
    return(NULL)
  }

  locus <- parts[1]
  fields <- parts[2]

  # Extract two-field (first two colon-separated fields)
  field_parts <- strsplit(fields, ":")[[1]]
  if (length(field_parts) < 1) {
    return(NULL)
  }

  # Handle P-group notation (e.g., "01:01P")
  if (length(field_parts) >= 2) {
    allele_2f <- paste(field_parts[1:2], collapse = ":")
  } else {
    # Single field (e.g., just "01")
    allele_2f <- field_parts[1]
  }

  list(locus = locus, allele_2f = allele_2f)
}

#' @noRd
.getSerologyPrefix <- function(locus) {
  # Map HLA loci to their serological prefixes
  prefix_map <- c(
    "A" = "A",
    "B" = "B",
    "C" = "Cw",
    "DRA" = "DR",
    "DRB1" = "DR",
    "DRB3" = "DR",
    "DRB4" = "DR",
    "DRB5" = "DR",
    "DQA1" = "DQ",
    "DQB1" = "DQ",
    "DPA1" = "DP",
    "DPB1" = "DP"
  )

  if (locus %in% names(prefix_map)) {
    return(prefix_map[locus])
  }

  # Default: use locus as-is
  locus
}

#' @noRd
.lookupSerology <- function(locus, allele_2f_input, serology_db) {
  # WMDA data has locus with asterisk (e.g., "A*", "DRB1*")
  locus_key <- paste0(locus, "*")

  # Direct lookup using data.table's native column reference in i
  # Need to avoid variable name conflicts with column names
  matches <- serology_db[locus == locus_key & allele_2f == allele_2f_input, ]

  if (nrow(matches) > 0 && !is.na(matches[["serology"]][1])) {
    return(matches[["serology"]][1])
  }

  # Try without leading zeros (e.g., "1:01" instead of "01:01")
  allele_2f_no_lead <- sub("^0+", "", allele_2f_input)
  matches <- serology_db[locus == locus_key & allele_2f == allele_2f_no_lead, ]

  if (nrow(matches) == 0) return(NA_character_)
  matches[["serology"]][1]
}

#' @noRd
.resolvePGroup <- function(locus_input, allele_2f_input, pgroups_db) {
  # P-groups data has locus WITHOUT asterisk (e.g., "A" not "A*")
  locus_key <- locus_input

  # Remove trailing P for lookup
  p_group_key <- sub("P$", "", allele_2f_input)

  # Use data.table's native column reference in i
  matches <- pgroups_db[locus == locus_key & p_group == p_group_key, ]

  if (nrow(matches) == 0 || is.na(matches[["reference_2f"]][1])) {
    # Try with P suffix in lookup
    matches <- pgroups_db[locus == locus_key & p_group == allele_2f_input, ]
  }

  if (nrow(matches) == 0 || is.na(matches[["reference_2f"]][1])) {
    return(NULL)
  }

  matches[["reference_2f"]][1]
}

#' @noRd
.resolveBroadToSplit <- function(input_locus, allele_2f_input, serology_input, wmda_data) {
  # Get the serology locus prefix
  ser_locus <- .getSerologyPrefix(input_locus)

  # Check if this serology is a broad antigen
  splits_db <- wmda_data$splits
  # Use data.table's native column reference in i
  matches <- splits_db[locus == ser_locus & broad == serology_input, ]

  if (nrow(matches) == 0 || is.na(matches[["splits"]][1])) {
    return(serology_input)
  }

  # Parse the available splits
  available_splits <- strsplit(matches[["splits"]][1], "\\|")[[1]]

  # Try to determine which split based on allele
  # Look up each split's allele patterns in serology database
  serology_db <- wmda_data$serology
  locus_key <- paste0(input_locus, "*")

  for (split in available_splits) {
    # Find alleles that map to this split
    split_matches <- serology_db[locus == locus_key & serology == split, ]
    split_alleles <- split_matches[["allele_2f"]]

    if (allele_2f_input %in% split_alleles) {
      return(split)
    }
  }

  # Could not determine split, return original serology
  serology_input
}

#' @noRd
.handleUnknown <- function(allele, na_action, message) {
  if (na_action == "error") {
    stop(message)
  }

  if (na_action == "warn") {
    warning(message)
  }

  data.table::data.table(
    allele = allele,
    locus = NA_character_,
    allele_2f = NA_character_,
    serology = NA_character_,
    serology_full = NA_character_
  )
}

#' @noRd
.loadWmdaData <- function() {
  # Check for cached/updated data first
  cache_dir <- .getWmdaCacheDir()
  cache_file <- file.path(cache_dir, "wmda_cache.rds")

  if (file.exists(cache_file)) {
    tryCatch({
      cached <- readRDS(cache_file)
      # Ensure keys are set (may be lost after loading)
      data.table::setkey(cached$serology, locus, allele_2f)
      data.table::setkey(cached$splits, locus, broad)
      data.table::setkey(cached$pgroups, locus, p_group)
      return(cached)
    }, error = function(e) {
      # Cache corrupted, fall through to bundled data
    })
  }

  # Access bundled package data - try multiple methods
  serology <- NULL
  splits <- NULL
  pgroups <- NULL

  # Method 1: Try from package namespace (when package is loaded)
  tryCatch({
    ns <- asNamespace("deepMatchR")
    serology <- get0("deepMatchR_wmda_serology", envir = ns, ifnotfound = NULL)
    splits <- get0("deepMatchR_wmda_splits", envir = ns, ifnotfound = NULL)
    pgroups <- get0("deepMatchR_wmda_pgroups", envir = ns, ifnotfound = NULL)
  }, error = function(e) NULL)

  # Method 2: Try using utils::data() to load from package
  if (is.null(serology) || is.null(splits) || is.null(pgroups)) {
    tryCatch({
      env <- new.env()
      data("deepMatchR_wmda_serology", package = "deepMatchR", envir = env)
      data("deepMatchR_wmda_splits", package = "deepMatchR", envir = env)
      data("deepMatchR_wmda_pgroups", package = "deepMatchR", envir = env)
      serology <- env$deepMatchR_wmda_serology
      splits <- env$deepMatchR_wmda_splits
      pgroups <- env$deepMatchR_wmda_pgroups
    }, error = function(e) NULL)
  }

  if (is.null(serology) || is.null(splits) || is.null(pgroups)) {
    stop("WMDA data not found. Package data may not be properly installed.")
  }

  # Make copies and set keys (keys are lost during lazy loading)
  serology <- data.table::copy(serology)
  splits <- data.table::copy(splits)
  pgroups <- data.table::copy(pgroups)

  data.table::setkey(serology, locus, allele_2f)
  data.table::setkey(splits, locus, broad)
  data.table::setkey(pgroups, locus, p_group)

  list(
    serology = serology,
    splits = splits,
    pgroups = pgroups
  )
}

#' @noRd
.getWmdaCacheDir <- function() {
  cache_base <- Sys.getenv("DEEPMATCHR_CACHE_DIR", "")
  if (nzchar(cache_base)) {
    return(cache_base)
  }

  # Default: use R user data directory
  rappdirs_available <- requireNamespace("rappdirs", quietly = TRUE)
  if (rappdirs_available) {
    return(rappdirs::user_cache_dir("deepMatchR"))
  }

  # Fallback
  file.path(tempdir(), "deepMatchR_cache")
}

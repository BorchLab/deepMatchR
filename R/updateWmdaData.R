#' Update WMDA Nomenclature Data
#'
#' Downloads fresh WMDA nomenclature data from the IMGT/HLA GitHub repository
#' and caches it locally for use by \code{\link{toSerology}}.
#'
#' @param version Character. The IMGT/HLA version to download. Default is "Latest"
#'   for the most recent release. Can also be a specific version tag (e.g., "3.54.0").
#' @param cache_dir Character. Directory to store cached data. If NULL (default),
#'   uses the package's default cache location.
#' @param force Logical. If TRUE, re-download even if cached data exists.
#' @param verbose Logical. If TRUE (default), print progress messages.
#'
#' @return Invisibly returns the path to the cache directory.
#'
#' @details
#' The WMDA nomenclature files are downloaded from:
#' \url{https://github.com/ANHIG/IMGTHLA}
#'
#' Downloaded files:
#' \itemize{
#'   \item \code{rel_dna_ser.txt}: DNA to serology mappings
#'   \item \code{rel_ser_ser.txt}: Broad to split relationships
#'   \item \code{hla_nom_p.txt}: P-group definitions
#' }
#'
#' The cache is stored as a single RDS file for fast loading. To clear the cache
#' and revert to bundled data, delete the cache directory or set the environment
#' variable \code{DEEPMATCHR_CACHE_DIR} to a new location.
#'
#' @examples
#' \dontrun{
#' # Update to latest WMDA data
#' updateWmdaData()
#'
#' # Update to a specific version
#' updateWmdaData(version = "3.54.0")
#'
#' # Force re-download
#' updateWmdaData(force = TRUE)
#' }
#'
#' @seealso \code{\link{toSerology}}
#' @export
updateWmdaData <- function(version = "Latest",
                           cache_dir = NULL,
                           force = FALSE,
                           verbose = TRUE) {

  # Determine cache directory
  if (is.null(cache_dir)) {
    cache_dir <- .getWmdaCacheDir()
  }

  # Create cache directory if needed
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  cache_file <- file.path(cache_dir, "wmda_cache.rds")
  version_file <- file.path(cache_dir, "wmda_version.txt")

  # Check if update needed
  if (!force && file.exists(cache_file) && file.exists(version_file)) {
    cached_version <- readLines(version_file, n = 1, warn = FALSE)
    if (cached_version == version) {
      if (verbose) message("Cache is up to date (version: ", version, ")")
      return(invisible(cache_dir))
    }
  }

  # Build URL base
  base_url <- paste0("https://raw.githubusercontent.com/ANHIG/IMGTHLA/", version, "/wmda/")

  if (verbose) message("Downloading WMDA data (version: ", version, ")...")

  # Download and parse files
  tryCatch({
    # Download rel_dna_ser.txt
    if (verbose) message("  Downloading rel_dna_ser.txt...")
    dna_ser_lines <- .downloadWmdaFile("rel_dna_ser.txt", base_url)
    serology_dt <- .parseDnaSer(dna_ser_lines)

    # Download rel_ser_ser.txt
    if (verbose) message("  Downloading rel_ser_ser.txt...")
    ser_ser_lines <- .downloadWmdaFile("rel_ser_ser.txt", base_url)
    splits_dt <- .parseSerSer(ser_ser_lines)

    # Download hla_nom_p.txt
    if (verbose) message("  Downloading hla_nom_p.txt...")
    nom_p_lines <- .downloadWmdaFile("hla_nom_p.txt", base_url)
    pgroups_dt <- .parseNomP(nom_p_lines)

    # Set keys for fast lookup
    data.table::setkey(serology_dt, locus, allele_2f)
    data.table::setkey(splits_dt, locus, broad)
    data.table::setkey(pgroups_dt, locus, p_group)

    # Save to cache
    if (verbose) message("  Saving to cache...")
    wmda_data <- list(
      serology = serology_dt,
      splits = splits_dt,
      pgroups = pgroups_dt
    )
    saveRDS(wmda_data, cache_file)
    writeLines(version, version_file)

    if (verbose) {
      message("Done! Cached ", nrow(serology_dt), " serology mappings, ",
              nrow(splits_dt), " split mappings, ",
              nrow(pgroups_dt), " P-groups")
      message("Cache location: ", cache_dir)
    }

  }, error = function(e) {
    stop("Failed to update WMDA data: ", e$message)
  })

  invisible(cache_dir)
}

#' @noRd
.downloadWmdaFile <- function(filename, base_url) {
  url <- paste0(base_url, filename)

  lines <- tryCatch(
    readLines(url, warn = FALSE),
    error = function(e) {
      stop("Failed to download ", filename, " from ", url, ": ", e$message)
    }
  )

  # Remove comment lines (starting with #)
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nzchar(trimws(lines))]

  lines
}

#' @noRd
.parseDnaSer <- function(lines) {
  # Split by semicolon
  parsed <- strsplit(lines, ";")

  dt <- data.table::data.table(
    locus = vapply(parsed, `[`, character(1), 1),
    allele = vapply(parsed, `[`, character(1), 2),
    unambiguous = vapply(parsed, function(x) if (length(x) >= 3) x[3] else "", character(1)),
    possible = vapply(parsed, function(x) if (length(x) >= 4) x[4] else "", character(1)),
    assumed = vapply(parsed, function(x) if (length(x) >= 5) x[5] else "", character(1)),
    expert = vapply(parsed, function(x) if (length(x) >= 6) x[6] else "", character(1))
  )

  # Clean up empty strings to NA
  dt[unambiguous == "", unambiguous := NA_character_]
  dt[possible == "", possible := NA_character_]
  dt[assumed == "", assumed := NA_character_]
  dt[expert == "", expert := NA_character_]

  # Extract two-field allele
  dt[, allele_2f := sub("^([^:]+:[^:]+).*", "\\1", allele)]

  # Choose best serology
  dt[, serology := data.table::fifelse(
    !is.na(unambiguous), unambiguous,
    data.table::fifelse(!is.na(possible), possible,
                        data.table::fifelse(!is.na(assumed), assumed, expert))
  )]

  # Handle multiple assignments (take first)
  dt[, serology := sub("/.*", "", serology)]

  # Keep only rows with serology
  dt <- dt[!is.na(serology) & nzchar(serology)]

  # Deduplicate
  dt <- unique(dt, by = c("locus", "allele_2f"))

  dt[, .(locus, allele_2f, serology)]
}

#' @noRd
.parseSerSer <- function(lines) {
  parsed <- strsplit(lines, ";")

  dt <- data.table::data.table(
    locus = vapply(parsed, `[`, character(1), 1),
    broad = vapply(parsed, function(x) if (length(x) >= 2) x[2] else "", character(1)),
    split = vapply(parsed, function(x) if (length(x) >= 3) x[3] else "", character(1))
  )

  dt <- dt[nzchar(broad) & nzchar(split)]

  # Aggregate splits
  dt_agg <- dt[, .(splits = paste(unique(split), collapse = "|")), by = .(locus, broad)]

  dt_agg
}

#' @noRd
.parseNomP <- function(lines) {
  parsed <- strsplit(lines, ";")

  results <- lapply(parsed, function(x) {
    if (length(x) < 2) return(NULL)

    allele_part <- x[1]
    p_group <- x[2]

    if (!grepl("\\*", allele_part)) return(NULL)

    locus <- sub("\\*.*", "", allele_part)
    alleles_str <- sub("^[^*]+\\*", "", allele_part)
    alleles <- strsplit(alleles_str, "/")[[1]]
    alleles_2f <- sub("^([^:]+:[^:]+).*", "\\1", alleles)

    data.table::data.table(
      locus = locus,
      p_group = p_group,
      reference_2f = alleles_2f[1]
    )
  })

  dt <- data.table::rbindlist(results[!vapply(results, is.null, logical(1))])
  unique(dt, by = c("locus", "p_group"))
}

#' Clear WMDA Cache
#'
#' Removes the cached WMDA data, causing \code{\link{toSerology}} to revert
#' to using the bundled package data.
#'
#' @param cache_dir Character. Directory containing cached data. If NULL,
#'   uses the default cache location.
#' @param verbose Logical. If TRUE (default), print status message.
#'
#' @return Invisibly returns TRUE if cache was cleared, FALSE if no cache existed.
#'
#' @examples
#' \dontrun{
#' # Clear cached WMDA data
#' clearWmdaCache()
#' }
#'
#' @seealso \code{\link{updateWmdaData}}, \code{\link{toSerology}}
#' @export
clearWmdaCache <- function(cache_dir = NULL, verbose = TRUE) {
  if (is.null(cache_dir)) {
    cache_dir <- .getWmdaCacheDir()
  }

  cache_file <- file.path(cache_dir, "wmda_cache.rds")
  version_file <- file.path(cache_dir, "wmda_version.txt")

  removed <- FALSE

  if (file.exists(cache_file)) {
    file.remove(cache_file)
    removed <- TRUE
  }

  if (file.exists(version_file)) {
    file.remove(version_file)
    removed <- TRUE
  }

  if (verbose) {
    if (removed) {
      message("WMDA cache cleared")
    } else {
      message("No WMDA cache found")
    }
  }

  invisible(removed)
}

# Build WMDA Serology Data Files
# This script downloads and processes WMDA nomenclature files to create
# the bundled serology lookup tables for deepMatchR.
#
# Run this script from the package root directory:
#   source("data-raw/build-wmda-serology.R")

library(data.table)

# --- Configuration ---
WMDA_BASE_URL <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/"

# --- Helper Functions ---

download_wmda_file <- function(filename, base_url = WMDA_BASE_URL) {
  url <- paste0(base_url, filename)
  message("Downloading: ", url)

  lines <- tryCatch(
    readLines(url, warn = FALSE),
    error = function(e) {
      stop("Failed to download ", filename, ": ", e$message)
    }
  )

  # Remove comment lines (starting with #)
  lines <- lines[!grepl("^#", lines)]
  lines <- lines[nzchar(trimws(lines))]

  lines
}

# --- Parse rel_dna_ser.txt ---
# Format: Locus;Allele;Unambiguous;Possible;Assumed;Expert
# We want the "Unambiguous" or best available serology

parse_dna_ser <- function(lines) {
  # Split by semicolon
  parsed <- strsplit(lines, ";")

  dt <- data.table(
    locus = sapply(parsed, `[`, 1),
    allele = sapply(parsed, `[`, 2),
    unambiguous = sapply(parsed, `[`, 3),
    possible = sapply(parsed, `[`, 4),
    assumed = sapply(parsed, `[`, 5),
    expert = sapply(parsed, `[`, 6)
  )

  # Clean up NA values
  dt[unambiguous == "", unambiguous := NA_character_]
  dt[possible == "", possible := NA_character_]
  dt[assumed == "", assumed := NA_character_]
  dt[expert == "", expert := NA_character_]

  # Extract two-field allele (first two fields)
  dt[, allele_2f := sub("^([^:]+:[^:]+).*", "\\1", allele)]

  # Choose best serology: unambiguous > possible > assumed > expert
  dt[, serology := fifelse(
    !is.na(unambiguous), unambiguous,
    fifelse(!is.na(possible), possible,
            fifelse(!is.na(assumed), assumed, expert))
  )]

  # Handle multiple serology assignments (take first if multiple separated by /)
  dt[, serology := sub("/.*", "", serology)]

  # Keep only rows with serology
  dt <- dt[!is.na(serology) & nzchar(serology)]

  # Deduplicate by locus + allele_2f, keeping first (most specific) entry
  dt <- unique(dt, by = c("locus", "allele_2f"))

  # Return simplified table
  dt[, .(locus, allele_2f, serology)]
}

# --- Parse rel_ser_ser.txt ---
# Format: Locus;Broad;Split
# Maps broad antigens to their splits

parse_ser_ser <- function(lines) {
  parsed <- strsplit(lines, ";")

  dt <- data.table(
    locus = sapply(parsed, `[`, 1),
    broad = sapply(parsed, `[`, 2),
    split = sapply(parsed, `[`, 3)
  )

  # Clean up
  dt <- dt[!is.na(broad) & nzchar(broad)]
  dt <- dt[!is.na(split) & nzchar(split)]

  # Aggregate splits for each broad
  dt_agg <- dt[, .(splits = paste(unique(split), collapse = "|")), by = .(locus, broad)]

  dt_agg
}

# --- Parse hla_nom_p.txt ---
# Format: Locus*Allele_list;P-group
# P-groups are alleles with identical protein sequences

parse_nom_p <- function(lines) {
  # Split by semicolon
  parsed <- strsplit(lines, ";")

  results <- lapply(parsed, function(x) {
    if (length(x) < 2) return(NULL)

    allele_part <- x[1]
    p_group <- x[2]

    # Extract locus from allele_part (format: "A*01:01/01:02/...")
    if (!grepl("\\*", allele_part)) return(NULL)

    locus <- sub("\\*.*", "", allele_part)
    alleles_str <- sub("^[^*]+\\*", "", allele_part)

    # Split by / to get individual alleles
    alleles <- strsplit(alleles_str, "/")[[1]]

    # Extract two-field for each
    alleles_2f <- sub("^([^:]+:[^:]+).*", "\\1", alleles)

    data.table(
      locus = locus,
      p_group = p_group,
      reference_2f = alleles_2f[1]  # First allele is the reference
    )
  })

  dt <- rbindlist(results[!sapply(results, is.null)])
  unique(dt, by = c("locus", "p_group"))
}

# --- Main Build Process ---

message("Building WMDA serology data files...")

# Download files
dna_ser_lines <- download_wmda_file("rel_dna_ser.txt")
ser_ser_lines <- download_wmda_file("rel_ser_ser.txt")
nom_p_lines <- download_wmda_file("hla_nom_p.txt")

# Parse files
message("Parsing rel_dna_ser.txt...")
deepMatchR_wmda_serology <- parse_dna_ser(dna_ser_lines)
message("  Found ", nrow(deepMatchR_wmda_serology), " DNA-to-serology mappings")

message("Parsing rel_ser_ser.txt...")
deepMatchR_wmda_splits <- parse_ser_ser(ser_ser_lines)
message("  Found ", nrow(deepMatchR_wmda_splits), " broad-to-split mappings")

message("Parsing hla_nom_p.txt...")
deepMatchR_wmda_pgroups <- parse_nom_p(nom_p_lines)
message("  Found ", nrow(deepMatchR_wmda_pgroups), " P-group definitions")

# Set keys for fast lookups
setkey(deepMatchR_wmda_serology, locus, allele_2f)
setkey(deepMatchR_wmda_splits, locus, broad)
setkey(deepMatchR_wmda_pgroups, locus, p_group)

# Save data files directly (without usethis)
data_dir <- "data"
if (!dir.exists(data_dir)) dir.create(data_dir)

message("Saving data files...")
save(deepMatchR_wmda_serology, file = file.path(data_dir, "deepMatchR_wmda_serology.rda"), compress = "xz")
save(deepMatchR_wmda_splits, file = file.path(data_dir, "deepMatchR_wmda_splits.rda"), compress = "xz")
save(deepMatchR_wmda_pgroups, file = file.path(data_dir, "deepMatchR_wmda_pgroups.rda"), compress = "xz")

message("Done! Created:")
message("  - data/deepMatchR_wmda_serology.rda")
message("  - data/deepMatchR_wmda_splits.rda")
message("  - data/deepMatchR_wmda_pgroups.rda")

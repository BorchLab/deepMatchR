# ---- Set Up ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(usethis)
})
source("./R/utils.R")

# Helper: project-local path for raw CSVs
pp <- function(...) file.path("data-raw", "raw", ...)

# ---- Inputs ------------------------------------------------------------------
# Example dataset is shipped in /data
example_file <- system.file("data", "deepMatchR_example.rda", package = "deepMatchR")
stopifnot(file.exists(example_file))
load(example_file)  # loads: deepMatchR_example
stopifnot(exists("deepMatchR_example"))

# External registry CSVs (keep in version control under data-raw/raw/)
c1_file <- pp("EpletRegistry_ClassI.csv")
c2_file <- pp("EpletRegistry_ClassII.csv")
stopifnot(file.exists(c1_file), file.exists(c2_file))

# ---- Process example SAB/PRA -------------------------------------------------
SAB1 <- .processSAB(deepMatchR_example[[1]])
SAB2 <- .processSAB(deepMatchR_example[[2]])
PRA  <- .processPRA(deepMatchR_example[[3]])

all_alleles <- bind_rows(
  select(SAB1, allele, serology = antigen),
  select(SAB2, allele, serology = antigen),
  select(PRA,  allele, serology = antigen)
) %>%
  filter(!allele %in% c("Bw4", "Bw6")) %>%
  distinct()

# ---- Eplet dictionary --------------------------------------------------------
read_eplet_csv <- function(path) read_csv(path, show_col_types = FALSE)

ClassI <- read_eplet_csv(c1_file)[,-1]
ClassII <- read_eplet_csv(c2_file)[,-1]

to_long_eplet <- function(df) {
  df %>%
    separate_rows(Alleles, sep = ",\\s*") %>%
    mutate(Alleles = str_trim(Alleles)) %>%
    select(-description) %>%
    rename(
      eplet      = `Epitope Name`,
      exposition = exposition,
      reactivity = `Antibody Reactivity`,
      evidence   = evidence,
      allele     = Alleles
    ) %>%
    filter(!is.na(allele), allele != "")
}

Eplet_Dictionary <- bind_rows(to_long_eplet(ClassI), to_long_eplet(ClassII))

deepMatchR_eplets <- Eplet_Dictionary 

# ---- CREG mapping ------------------------------------------------------------
creg_tbl <- tibble::tribble(
  ~CREG, ~Serology,
  "1C",  "A1, 3, 9 (23, 24), 11, 29, 30, 31, 36, 80",
  "10C", "A10 (25, 26, 34, 66), 11, 28 (68, 69), 32, 33, 43, 74",
  "2C",  "A2, 9 (23, 24), 28 (68, 69), B17 (57, 58)",
  "5C",  "B5 (51, 52), 15 (62, 63, 75, 76, 77), 17 (57, 58), 18, 21 (49, 50), 35, 46, 53, 70 (71, 72), 73, 78",
  "7C",  "B7, 8, 13, 22 (54, 55, 56), 27, 40 (60, 61), 41, 42, 47, 48, 59, 67, 81, 82",
  "8C",  "B8, 14 (64, 65), 16 (38, 39), 18, 59, 67",
  "12C", "B12 (44, 45), 13, 21 (49, 50), 37, 40 (60, 61), 41, 47",
  "Bw4", "A23, 24, 25, 32, B13, 27, 37, 38, 44, 47, 49, 51, 52, 53, 57, 58, 59, 63, 77",
  "Bw6", "B7, 8, 18, 35, 39, 41, 42, 45, 46, 48, 50, 54, 55, 56, 60, 61, 62, 64, 65, 67, 71, 72, 73, 75, 76, 78, 81, 82"
)

# Split on commas not inside parentheses
split_top_level_commas <- function(s) {
  out <- character(); buf <- ""; depth <- 0L
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  for (ch in chars) {
    if (ch == "(") { depth <- depth + 1L; buf <- paste0(buf, ch)
    } else if (ch == ")") { depth <- depth - 1L; buf <- paste0(buf, ch)
    } else if (ch == "," && depth == 0L) { out <- c(out, buf); buf <- ""
    } else { buf <- paste0(buf, ch) }
  }
  c(out, buf) |> str_trim() |> discard(~ .x == "")
}

expand_serology <- function(s) {
  items <- split_top_level_commas(s)
  out <- character(); prefix <- NA_character_
  for (tok in items) {
    tok <- str_squish(tok)
    m <- str_match(tok, "^(A|B)(\\d+)\\s*\\(([^)]*)\\)$")
    if (!anyNA(m[1, ])) {
      prefix <- m[1, 2]; base_num <- m[1, 3]
      inner <- str_split(m[1, 4], ",", simplify = TRUE) |> str_trim()
      inner <- inner[nchar(inner) > 0]
      out <- c(out, paste0(prefix, base_num), paste0(prefix, inner)); next
    }
    m <- str_match(tok, "^(\\d+)\\s*\\(([^)]*)\\)$")
    if (!anyNA(m[1, ])) {
      if (is.na(prefix)) warning("No prefix in context for token: ", tok)
      base_num <- m[1, 2]
      inner <- str_split(m[1, 3], ",", simplify = TRUE) |> str_trim()
      inner <- inner[nchar(inner) > 0]
      out <- c(out, paste0(prefix, base_num), paste0(prefix, inner)); next
    }
    m <- str_match(tok, "^(A|B)(\\d+)$")
    if (!anyNA(m[1, ])) { prefix <- m[1, 2]; out <- c(out, paste0(prefix, m[1, 3])); next }
    m <- str_match(tok, "^(\\d+)$")
    if (!anyNA(m[1, ])) { if (is.na(prefix)) warning("No prefix for: ", tok); out <- c(out, paste0(prefix, m[1, 2])); next }
    warning("Unrecognized token: ", tok)
  }
  unique(out)
}

creg_map <- creg_tbl %>%
  mutate(serology_list = map(Serology, expand_serology)) %>%
  select(CREG, serology_list) %>%
  unnest_longer(serology_list, values_to = "serology") %>%
  distinct(serology, CREG)

# Join allele↔serology with serology↔CREG to get allele↔CREG
deepMatchR_cregs <- all_alleles %>%
  left_join(creg_map, by = "serology", relationship = "many-to-many") %>%
  arrange(allele, CREG) %>%
  distinct()

# ---- Save data (compressed) ---------------------------------------------------
usethis::use_data(deepMatchR_eplets, deepMatchR_cregs, overwrite = TRUE, compress = "xz")



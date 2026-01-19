#' Calculate Peptide Binding Load for Transplant Risk Assessment
#'
#' @description
#' Predicts transplant risk by calculating peptide-HLA binding affinities
#' between recipient HLA molecules and donor-mismatched peptides. Supports
#' multiple binding prediction backends: built-in position weight matrix (PWM),
#' NetMHCpan, or MHCflurry.
#'
#' @param recipient An `hla_genotype` object or character vector of HLA allele names.
#' @param donor An `hla_genotype` object, character vector of HLA allele names, or
#'   a character vector of peptide sequences. If `hla_genotype` or allele names,
#'   mismatched peptides are derived automatically from sequence differences.
#' @param backend Character. Binding prediction method: `"pwm"` (default, no external
#'   dependencies), `"netmhcpan"`, or `"mhcflurry"`.
#' @param backend_path Character. Path to external tool executable. Required for
#'   `"netmhcpan"` backend.
#' @param peptide_length Integer. Peptide length(s) to consider. Default `9L`.
#' @param binding_threshold Numeric. IC50 threshold (nM) for "strong binder".
#'   Default `500`.
#' @param weak_threshold Numeric. IC50 threshold (nM) for "weak binder".
#'   Default `5000`.
#' @param return Character. One of `"total"` (risk score), `"summary"` (per-allele),
#'   or `"detail"` (per-peptide table). Default `"total"`.
#' @param aggregate_method Character. How to combine per-peptide scores:
#'   `"sum"`, `"max"`, `"mean"`. Default `"sum"`.
#'
#' @return Depends on `return`:
#'   - `"total"`: Numeric risk score.
#'   - `"summary"`: data.frame with per-HLA-allele binding summary.
#'   - `"detail"`: data.frame with columns: `peptide`, `hla_allele`, `predicted_ic50`,
#'     `binding_level`, `contribution`.
#'
#' @details
#' The function works in several steps:
#' 1. **Input processing**: Converts inputs to standard format (allele names and peptides)
#' 2. **Peptide derivation**: If donor is genotype/alleles, derives mismatched peptides
#'    by comparing sequences and generating overlapping k-mers from mismatch regions

#' 3. **Binding prediction**: Uses selected backend to predict IC50 values
#' 4. **Risk calculation**: Aggregates binding predictions into a risk score
#'
#' The **PWM backend** uses simplified position weight matrices based on HLA supertypes.
#' For production use with high accuracy requirements, NetMHCpan is recommended.
#'
#' Risk score formula:
#' \deqn{contribution = (1 - IC50/weak\_threshold) \times multiplier}
#' where multiplier is 2 for strong binders, 1 for weak binders.
#'
#' @examples
#' # Create donor and recipient genotypes
#' recipient <- data.frame(
#'   A_1 = "A*02:01", A_2 = "A*03:01",
#'   B_1 = "B*07:02", B_2 = "B*44:02"
#' )
#' donor <- data.frame(
#'   A_1 = "A*01:01", A_2 = "A*24:02",
#'   B_1 = "B*08:01", B_2 = "B*35:01"
#' )
#' rgeno <- hlaGeno(recipient)
#' dgeno <- hlaGeno(donor)
#'
#' \dontrun{
#' # Calculate total binding load (requires sequence data)
#' calculatePeptideBindingLoad(rgeno, dgeno)
#'
#' # Get detailed per-peptide results
#' calculatePeptideBindingLoad(rgeno, dgeno, return = "detail")
#'
#' # Use with raw peptides
#' peptides <- c("GILGFVFTL", "NLVPMVATV", "FLKEKGGL")
#' calculatePeptideBindingLoad(rgeno, peptides)
#' }
#'
#' @seealso \code{\link{calculateMismatchLoad}}, \code{\link{quantifyMismatch}}
#'
#' @export
calculatePeptideBindingLoad <- function(
    recipient,
    donor,
    backend = c("pwm", "netmhcpan", "mhcflurry"),
    backend_path = NULL,
    peptide_length = 9L,
    binding_threshold = 500,
    weak_threshold = 5000,
    return = c("total", "summary", "detail"),
    aggregate_method = c("sum", "max", "mean")
) {
  backend <- match.arg(backend)
  return <- match.arg(return)
  aggregate_method <- match.arg(aggregate_method)

 # --- 1. Process recipient to get HLA alleles ---
  recipient_alleles <- .extractAlleles(recipient)
  if (length(recipient_alleles) == 0L) {
    stop("No valid recipient HLA alleles found.")
  }

  # --- 2. Process donor to get peptides ---
  peptides <- .getPeptides(donor, recipient, peptide_length)
  if (length(peptides) == 0L) {
    if (return == "total") return(0)
    if (return == "summary") {
      return(data.frame(
        hla_allele = recipient_alleles,
        n_peptides = 0L,
        n_strong = 0L,
        n_weak = 0L,
        risk_contribution = 0,
        stringsAsFactors = FALSE
      ))
    }
    return(data.frame(
      peptide = character(0),
      hla_allele = character(0),
      predicted_ic50 = numeric(0),
      binding_level = character(0),
      contribution = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  # --- 3. Predict binding for each peptide x allele combination ---
  results <- .predictBinding(
    peptides = peptides,
    alleles = recipient_alleles,
    backend = backend,
    backend_path = backend_path,
    peptide_length = peptide_length
  )

  # --- 4. Calculate binding levels and contributions ---
  results$binding_level <- ifelse(
    results$predicted_ic50 <= binding_threshold, "strong",
    ifelse(results$predicted_ic50 <= weak_threshold, "weak", "non_binder")
  )

  results$contribution <- ifelse(
    results$binding_level == "strong",
    2 * (1 - results$predicted_ic50 / weak_threshold),
    ifelse(
      results$binding_level == "weak",
      1 * (1 - results$predicted_ic50 / weak_threshold),
      0
    )
  )
  results$contribution <- pmax(results$contribution, 0)

  # --- 5. Return based on requested format ---
  if (return == "detail") {
    return(results)
  }

  if (return == "summary") {
    summary_df <- do.call(rbind, lapply(recipient_alleles, function(a) {
      subset <- results[results$hla_allele == a, ]
      data.frame(
        hla_allele = a,
        n_peptides = nrow(subset),
        n_strong = sum(subset$binding_level == "strong"),
        n_weak = sum(subset$binding_level == "weak"),
        risk_contribution = sum(subset$contribution),
        stringsAsFactors = FALSE
      )
    }))
    return(summary_df)
  }

  # return == "total"
  total_risk <- switch(aggregate_method,
    sum = sum(results$contribution),
    max = if (nrow(results) > 0) max(results$contribution) else 0,
    mean = if (nrow(results) > 0) mean(results$contribution) else 0
  )
  return(total_risk)
}


# --- Internal helper functions ---

#' Extract allele names from various input types
#' @keywords internal
.extractAlleles <- function(x) {
  if (inherits(x, "hla_genotype")) {
    validateHlaGeno(x)
    alleles <- unlist(x$data[1, , drop = TRUE])
    alleles <- alleles[!is.na(alleles) & nzchar(alleles)]
    return(unique(as.character(alleles)))
  }
  if (is.character(x)) {
    # Check if these look like alleles (contain *)
    if (all(grepl("\\*", x))) {
      return(unique(x))
    }
  }
  stop("Input must be an hla_genotype object or character vector of allele names.")
}


#' Get peptides from donor input
#' @keywords internal
.getPeptides <- function(donor, recipient, peptide_length) {
  # If donor is already peptides (character without *)
  if (is.character(donor) && !any(grepl("\\*", donor))) {
    # Filter to requested peptide length
    peptides <- donor[nchar(donor) == peptide_length]
    return(unique(peptides))
  }

  # Otherwise, derive peptides from sequence comparison
  donor_alleles <- .extractAlleles(donor)
  recipient_alleles <- .extractAlleles(recipient)

  # Get mismatched positions and generate peptides
  peptides <- .deriveMismatchedPeptides(
    donor_alleles = donor_alleles,
    recipient_alleles = recipient_alleles,
    peptide_length = peptide_length
  )
  return(unique(peptides))
}


#' Derive mismatched peptides from sequence comparison
#' @keywords internal
.deriveMismatchedPeptides <- function(donor_alleles, recipient_alleles, peptide_length) {
  peptides <- character(0)

  for (d_allele in donor_alleles) {
    # Get donor sequence
    d_seq <- tryCatch(
      getAlleleSequence(d_allele),
      error = function(e) NULL
    )
    if (is.null(d_seq) || nchar(d_seq) < peptide_length) next

    # Compare against each recipient allele of same locus
    d_locus <- sub("\\*.*", "", d_allele)

    for (r_allele in recipient_alleles) {
      r_locus <- sub("\\*.*", "", r_allele)
      if (d_locus != r_locus) next

      r_seq <- tryCatch(
        getAlleleSequence(r_allele),
        error = function(e) NULL
      )
      if (is.null(r_seq)) next

      # Find mismatch positions
      mismatch_detail <- quantifyMismatch(r_seq, d_seq, return = "detail")
      if (nrow(mismatch_detail) == 0) next

      # Generate peptides around mismatch positions
      mismatch_positions <- mismatch_detail$position

      for (pos in mismatch_positions) {
        # Generate all peptides that include this position
        start_min <- max(1, pos - peptide_length + 1)
        start_max <- min(pos, nchar(d_seq) - peptide_length + 1)

        if (start_max >= start_min) {
          for (start in start_min:start_max) {
            pep <- substr(d_seq, start, start + peptide_length - 1)
            if (nchar(pep) == peptide_length && !grepl("-|X", pep)) {
              peptides <- c(peptides, pep)
            }
          }
        }
      }
    }
  }

  return(unique(peptides))
}


#' Predict binding using selected backend
#' @keywords internal
.predictBinding <- function(peptides, alleles, backend, backend_path, peptide_length) {
  if (backend == "pwm") {
    return(.predictBindingPWM(peptides, alleles))
  } else if (backend == "netmhcpan") {
    return(.predictBindingNetMHCpan(peptides, alleles, backend_path))
  } else if (backend == "mhcflurry") {
    return(.predictBindingMHCflurry(peptides, alleles))
  }
  stop("Unknown backend: ", backend)
}


#' PWM-based binding prediction (built-in)
#' @keywords internal
.predictBindingPWM <- function(peptides, alleles) {
  # Simplified PWM-based prediction using HLA supertype scoring
  # This provides a reasonable approximation without external dependencies

  # HLA supertype mapping (simplified)
  supertypes <- list(
    A02 = c("A*02:01", "A*02:02", "A*02:03", "A*02:06", "A*02:07", "A*02:11", "A*68:02"),
    A03 = c("A*03:01", "A*11:01", "A*31:01", "A*33:01", "A*68:01"),
    A01 = c("A*01:01", "A*26:01", "A*32:01"),
    A24 = c("A*24:02", "A*23:01"),
    B07 = c("B*07:02", "B*35:01", "B*51:01", "B*53:01", "B*54:01", "B*55:01"),
    B44 = c("B*44:02", "B*44:03", "B*18:01", "B*37:01", "B*40:01", "B*41:01"),
    B27 = c("B*27:05", "B*14:01", "B*38:01", "B*39:01"),
    B58 = c("B*58:01", "B*57:01"),
    B08 = c("B*08:01")
  )

  # Preferred anchor residues for each supertype (position 2 and C-terminus)
  # These are simplified for the PWM approach
  anchor_p2 <- list(
    A02 = c("L", "M", "V", "I", "A", "T"),
    A03 = c("L", "V", "M", "I", "S", "A", "T"),
    A01 = c("T", "S", "M", "L"),
    A24 = c("Y", "F", "W"),
    B07 = c("P"),
    B44 = c("E", "D"),
    B27 = c("R", "K", "H"),
    B58 = c("A", "S", "T"),
    B08 = c("K", "R")
  )

  anchor_c <- list(
    A02 = c("L", "V", "I", "M", "A"),
    A03 = c("K", "R", "Y"),
    A01 = c("Y"),
    A24 = c("F", "L", "I", "W"),
    B07 = c("L", "M", "F"),
    B44 = c("Y", "F", "W"),
    B27 = c("L", "F", "K", "R"),
    B58 = c("W", "F", "Y"),
    B08 = c("L")
  )

  results <- data.frame(
    peptide = character(0),
    hla_allele = character(0),
    predicted_ic50 = numeric(0),
    stringsAsFactors = FALSE
  )

  for (allele in alleles) {
    # Find supertype for this allele
    supertype <- NULL
    for (st in names(supertypes)) {
      # Check exact match or prefix match
      if (allele %in% supertypes[[st]] ||
          any(startsWith(allele, sub("\\*.*", "*", supertypes[[st]])))) {
        supertype <- st
        break
      }
    }

    # Default to A02-like if unknown
    if (is.null(supertype)) supertype <- "A02"

    for (pep in peptides) {
      if (nchar(pep) < 2) next

      # Score based on anchor positions
      p2 <- substr(pep, 2, 2)
      pc <- substr(pep, nchar(pep), nchar(pep))

      score <- 0

      # P2 anchor contribution
      if (p2 %in% anchor_p2[[supertype]]) {
        score <- score + 2
      } else if (p2 %in% c("L", "M", "V", "I", "A")) {
        score <- score + 1
      }

      # C-terminal anchor contribution
      if (pc %in% anchor_c[[supertype]]) {
        score <- score + 2
      } else if (pc %in% c("L", "V", "I", "F", "Y")) {
        score <- score + 1
      }

      # Convert score to approximate IC50
      # Higher score = lower IC50 (better binding)
      # Score 4 -> ~100nM, Score 0 -> ~10000nM
      ic50 <- 10000 / (2^score)
      ic50 <- max(50, min(50000, ic50))

      results <- rbind(results, data.frame(
        peptide = pep,
        hla_allele = allele,
        predicted_ic50 = ic50,
        stringsAsFactors = FALSE
      ))
    }
  }

  return(results)
}


#' NetMHCpan-based binding prediction
#' @keywords internal
.predictBindingNetMHCpan <- function(peptides, alleles, backend_path) {
  if (is.null(backend_path) || !file.exists(backend_path)) {
    stop("NetMHCpan backend requires 'backend_path' to point to the netMHCpan executable.")
  }

  results <- data.frame(
    peptide = character(0),
    hla_allele = character(0),
    predicted_ic50 = numeric(0),
    stringsAsFactors = FALSE
  )

  # Create temp file for peptides
  pep_file <- tempfile(fileext = ".pep")
  writeLines(peptides, pep_file)
  on.exit(unlink(pep_file), add = TRUE)

  for (allele in alleles) {
    # Format allele for NetMHCpan (e.g., HLA-A02:01)
    formatted_allele <- gsub("\\*", "", allele)
    formatted_allele <- paste0("HLA-", formatted_allele)

    # Run NetMHCpan
    out_file <- tempfile(fileext = ".out")
    on.exit(unlink(out_file), add = TRUE)

    cmd <- sprintf(
      "%s -p %s -a %s -BA > %s 2>&1",
      backend_path, pep_file, formatted_allele, out_file
    )

    system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)

    # Parse output
    if (file.exists(out_file)) {
      lines <- readLines(out_file, warn = FALSE)
      # Find data lines (they start with position numbers after header)
      data_lines <- grep("^\\s*[0-9]", lines, value = TRUE)

      for (line in data_lines) {
        fields <- strsplit(trimws(line), "\\s+")[[1]]
        if (length(fields) >= 12) {
          pep <- fields[3]
          ic50 <- as.numeric(fields[12])
          if (!is.na(ic50) && pep %in% peptides) {
            results <- rbind(results, data.frame(
              peptide = pep,
              hla_allele = allele,
              predicted_ic50 = ic50,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }

  # If NetMHCpan didn't return results for some peptides, add them with NA
  missing <- setdiff(
    paste(rep(peptides, each = length(alleles)), rep(alleles, length(peptides)), sep = "_"),
    paste(results$peptide, results$hla_allele, sep = "_")
  )

  if (length(missing) > 0) {
    missing_parts <- strsplit(missing, "_")
    for (mp in missing_parts) {
      results <- rbind(results, data.frame(
        peptide = mp[1],
        hla_allele = mp[2],
        predicted_ic50 = 50000,  # Assume non-binder
        stringsAsFactors = FALSE
      ))
    }
  }

  return(results)
}


#' MHCflurry-based binding prediction
#' @keywords internal
.predictBindingMHCflurry <- function(peptides, alleles) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for MHCflurry backend. ",
         "Install it with: install.packages('reticulate')")
  }

  # Check if mhcflurry is available
  mhcflurry_available <- tryCatch({
    reticulate::py_module_available("mhcflurry")
  }, error = function(e) FALSE)

  if (!mhcflurry_available) {
    stop("Python package 'mhcflurry' is not available. ",
         "Install it with: pip install mhcflurry && mhcflurry-downloads fetch")
  }

  # Import mhcflurry
  mhcflurry <- reticulate::import("mhcflurry")
  predictor <- mhcflurry$Class1PresentationPredictor$load()

  results <- data.frame(
    peptide = character(0),
    hla_allele = character(0),
    predicted_ic50 = numeric(0),
    stringsAsFactors = FALSE
  )

  # Format alleles for mhcflurry (e.g., HLA-A*02:01)
  formatted_alleles <- paste0("HLA-", alleles)

  # Run predictions
  for (i in seq_along(formatted_alleles)) {
    allele <- formatted_alleles[i]
    original_allele <- alleles[i]

    tryCatch({
      predictions <- predictor$predict(
        peptides = peptides,
        alleles = rep(allele, length(peptides))
      )

      pred_df <- reticulate::py_to_r(predictions)

      results <- rbind(results, data.frame(
        peptide = pred_df$peptide,
        hla_allele = original_allele,
        predicted_ic50 = pred_df$mhcflurry_affinity,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      # If prediction fails for this allele, add entries with high IC50
      results <- rbind(results, data.frame(
        peptide = peptides,
        hla_allele = original_allele,
        predicted_ic50 = 50000,
        stringsAsFactors = FALSE
      ))
    })
  }

  return(results)
}

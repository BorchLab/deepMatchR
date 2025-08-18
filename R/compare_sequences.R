#' Compare two HLA protein sequences and identify polymorphisms
#'
#' This function compares two HLA protein sequences of the same length and
#' identifies the positions at which the amino acids differ. It also
#' characterizes the change in physicochemical properties (polarity and charge)
#' at each polymorphic site based on IMGT classifications.
#'
#' @param seq1 A character string representing the first protein sequence.
#' @param seq2 A character string representing the second protein sequence.
#' @return A data frame with the following columns for each polymorphism:
#'   \itemize{
#'     \item \code{Position}: The position of the polymorphism (1-based index).
#'     \item \code{AA_Seq1}: The amino acid in sequence 1 at the polymorphic site.
#'     \item \code{AA_Seq2}: The amino acid in sequence 2 at the polymorphic site.
#'     \item \code{Polarity_Change}: A string describing the change in polarity (e.g., "Polar to Nonpolar").
#'     \item \code{Charge_Change}: A string describing the change in charge (e.g., "Positive to Negative").
#'   }
#' If the sequences are identical, it returns an empty data frame.
#' @examples
#' seq1 <- "YFAMYGEKVAHTHVDTLYVRYHY"
#' seq2 <- "YFDMYGEKVAHTHVDTLYVRYHY"
#' compare_hla_sequences(seq1, seq2)
#' @export
compare_hla_sequences <- function(seq1, seq2) {
  # Ensure sequences are character strings
  if (!is.character(seq1) || !is.character(seq2)) {
    stop("Input sequences must be character strings.")
  }

  # Ensure sequences are of the same length
  if (nchar(seq1) != nchar(seq2)) {
    stop("Input sequences must be of the same length.")
  }

  # Split sequences into character vectors
  s1 <- strsplit(seq1, "")[[1]]
  s2 <- strsplit(seq2, "")[[1]]

  # Define amino acid classifications based on IMGT
  polarity_map <- c(
    'R'='Polar', 'N'='Polar', 'D'='Polar', 'Q'='Polar', 'E'='Polar',
    'H'='Polar', 'K'='Polar', 'S'='Polar', 'T'='Polar', 'Y'='Polar',
    'A'='Nonpolar', 'C'='Nonpolar', 'G'='Nonpolar', 'I'='Nonpolar',
    'L'='Nonpolar', 'M'='Nonpolar', 'F'='Nonpolar', 'P'='Nonpolar',
    'W'='Nonpolar', 'V'='Nonpolar'
  )

  charge_map <- c(
    'R'='Positive', 'H'='Positive', 'K'='Positive',
    'D'='Negative', 'E'='Negative',
    'A'='Uncharged', 'N'='Uncharged', 'C'='Uncharged', 'Q'='Uncharged',
    'G'='Uncharged', 'I'='Uncharged', 'L'='Uncharged', 'M'='Uncharged',
    'F'='Uncharged', 'P'='Uncharged', 'S'='Uncharged', 'T'='Uncharged',
    'W'='Uncharged', 'Y'='Uncharged', 'V'='Uncharged'
  )

  # Find differing positions
  diff_indices <- which(s1 != s2)

  # If no differences, return an empty data frame
  if (length(diff_indices) == 0) {
    return(data.frame(
      Position = integer(),
      AA_Seq1 = character(),
      AA_Seq2 = character(),
      Polarity_Change = character(),
      Charge_Change = character(),
      stringsAsFactors = FALSE
    ))
  }

  # Function to get property change description
  get_change_desc <- function(val1, val2, property_map) {
    prop1 <- property_map[val1]
    prop2 <- property_map[val2]
    if (is.na(prop1) || is.na(prop2)) return("Unknown")
    if (prop1 == prop2) return("No change")
    return(paste(prop1, "to", prop2))
  }

  # Collect polymorphism data
  polymorphisms <- lapply(diff_indices, function(i) {
    aa1 <- s1[i]
    aa2 <- s2[i]

    # Check for unknown amino acids
    if (!aa1 %in% names(polarity_map) || !aa2 %in% names(polarity_map)) {
        warning(paste0("Unknown amino acid at position ", i, ". Skipping property analysis for this position."))
        polarity_change <- "Unknown"
        charge_change <- "Unknown"
    } else {
        polarity_change <- get_change_desc(aa1, aa2, polarity_map)
        charge_change <- get_change_desc(aa1, aa2, charge_map)
    }


    data.frame(
      Position = i,
      AA_Seq1 = aa1,
      AA_Seq2 = aa2,
      Polarity_Change = polarity_change,
      Charge_Change = charge_change,
      stringsAsFactors = FALSE
    )
  })

  # Combine list of data frames into one
  do.call(rbind, polymorphisms)
}

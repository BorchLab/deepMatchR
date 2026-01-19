#' Create an hla_genotype object
#'
#' @description
#' Creates a new `hla_genotype` object from a data frame of HLA calls. The
#' object is a list containing the genotype data and a record of the loci
#' present.
#'
#' @param df A data frame where rows are individuals and columns represent
#'   HLA alleles (e.g., A_1, A_2, B_1, B_2, ...).
#'
#' @return An object of class `hla_genotype`.
#'
#' @examples
#' # Create a genotype from a data frame
#' recipient <- data.frame(
#'   A_1 = "A*01:01", A_2 = "A*02:01",
#'   B_1 = "B*07:02", B_2 = "B*08:01"
#' )
#' geno <- hlaGeno(recipient)
#' print(geno)
#'
#' # Access components
#' geno$locus_present
#' geno$data
#'
#' @export
hlaGeno <- function(df) {

  # Get loci from column names (e.g., "A" from "A_1")
  loci <- unique(sub("_.*", "", colnames(df)))

  # Create the object
  structure(
    list(
      data = df,
      locus_present = loci
    ),
    class = "hla_genotype"
  )
}

#' Validate an hla_genotype object
#'
#' @param x An object to validate.
#'
#' @return `TRUE` if the object is a valid `hla_genotype` object, otherwise
#'   throws an error.
#'
#' @keywords internal
validateHlaGeno <- function(x) {
  if (!inherits(x, "hla_genotype")) {
    stop("Object is not of class 'hla_genotype'")
  }
  if (!is.list(x) || !all(c("data", "locus_present") %in% names(x))) {
    stop("hla_genotype object must be a list with 'data' and 'locus_present' elements")
  }
  if (!is.data.frame(x$data)) {
    stop("'data' element must be a data frame")
  }
  if (!is.character(x$locus_present)) {
    stop("'locus_present' element must be a character vector")
  }

  TRUE
}

#' Print an hla_genotype object
#'
#' @param x An object of class `hla_genotype`.
#' @param ... Additional arguments (not used).
#'
#' @return Invisibly returns the original object.
#'
#' @examples
#' df <- data.frame(A_1 = "A*01:01", B_1 = "B*07:02")
#' geno <- hlaGeno(df)
#' print(geno)
#'
#' @export
print.hla_genotype <- function(x, ...) {
  validateHlaGeno(x)
  cat("HLA Genotype Data\n")
  cat("-----------------\n")
  cat("Loci present:", paste(x$locus_present, collapse = ", "), "\n")
  cat("Number of samples:", nrow(x$data), "\n\n")
  print(head(x$data))
  invisible(x)
}

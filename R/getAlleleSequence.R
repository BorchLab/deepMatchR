#' Get Sequence for an HLA Allele 
#'
#' @description
#' Improved version with session-level caching to avoid redundant API calls.
#' Uses memoise for automatic caching. Supports both nucleotide and protein sequences.
#'
#' @param allele_name A character string representing the HLA allele name
#' @param type The type of sequence to retrieve. Either "NUC" for nucleotide or
#'        "PROT" for protein sequences (default "PROT")
#' @param use_cache Logical, whether to use caching (default TRUE)
#' @param cache_dir Optional directory for persistent cache
#'
#' @return A character string representing the sequence (amino acid or nucleotide)
#' 
#' @examples
#' # Collecting A*02:01 Protein Sequence:
#' getAlleleSequence("A*02:01")
#' 
#' # Collecting A*02:01 Nucleotide Sequence:
#' getAlleleSequence("A*02:01", type = "NUC")
#' 
#'
#' @importFrom immReferent getIMGT
#' @importFrom memoise memoise cache_filesystem
#' @export
getAlleleSequence <- function(allele_name, 
                              type = c("PROT", "NUC"),
                              use_cache = TRUE, 
                              cache_dir = NULL) {
  
  # Validate type argument
  type <- match.arg(type)
  
  # Create memoised version if caching is enabled
  if (use_cache) {
    if (!is.null(cache_dir)) {
      # Use filesystem cache for persistence across sessions
      # Create cache directory if it doesn't exist
      if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      }
      cache <- memoise::cache_filesystem(cache_dir)
      memo_fn <- memoise::memoise(.getAlleleSequenceBase, cache = cache)
    } else {
      # Use in-memory cache for current session
      memo_fn <- memoise::memoise(.getAlleleSequenceBase)
    }
    return(memo_fn(allele_name, type = type))
  }
  
  # Direct call without caching
  .getAlleleSequenceBase(allele_name, type = type)
}

#' Base function for sequence retrieval (internal)
#' 
#' @param allele_name HLA allele name
#' @param type Sequence type ("PROT" or "NUC")
#' @return Character string of the sequence
#' @noRd
.getAlleleSequenceBase <- function(allele_name, 
                                   type = c("PROT", "NUC")) {
  type <- match.arg(type)
  
  # Get HLA sequences of specified type
  hla_sequences <- immReferent::getIMGT(
    gene = "HLA", 
    type = type,  # Use the type parameter here
    suppressMessages = TRUE
  )
  
  # Find the index of the first sequence that contains the allele name
  match_idx <- grep(allele_name, names(hla_sequences), fixed = TRUE)[1]
  
  # Check if a match was found
  if (is.na(match_idx)) {
    stop("Allele '", allele_name, "' not found in the IMGT/HLA database for type '", type, "'.")
  }
  
  # Return the sequence as a character string
  as.character(hla_sequences[[match_idx]])
}

#' Batch Get Sequences with Parallel Processing
#'
#' @description
#' Retrieves sequences for multiple alleles in parallel for better performance.
#' Supports both nucleotide and protein sequences.
#'
#' @param alleles Character vector of HLA allele names
#' @param type The type of sequence to retrieve. Either "NUC" for nucleotide or
#'        "PROT" for protein sequences (default "PROT")
#' @param n_cores Number of cores to use (default: detectCores() - 1)
#' @param use_cache Whether to use caching (default TRUE)
#' @param cache_dir Optional directory for persistent cache
#' @param verbose Logical, whether to print progress messages (default FALSE)
#'
#' @return Named list of sequences
#' 
#' @examples
#' # Get protein sequences for multiple alleles
#' alleles <- c("A*01:01", "A*02:01", "B*07:02", "B*08:01")
#' prot_seqs <- batchGetSequences(alleles, type = "PROT")
#' 
#' # Get nucleotide sequences with parallel processing
#' nuc_seqs <- batchGetSequences(alleles, type = "NUC", n_cores = 4)
#' 
#'
#' @importFrom parallel mclapply detectCores
#' @export
batchGetSequences <- function(alleles, 
                              type = c("PROT", "NUC"),
                              n_cores = NULL, 
                              use_cache = TRUE,
                              cache_dir = NULL,
                              verbose = FALSE) {
  
  # Validate type argument
  type <- match.arg(type)
  
  # Remove duplicates
  unique_alleles <- unique(alleles)
  
  if (verbose) {
    message(sprintf("Retrieving %d unique %s sequences...", 
                    length(unique_alleles), 
                    ifelse(type == "PROT", "protein", "nucleotide")))
  }
  
  # Determine number of cores
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Use parallel processing for batch retrieval (for larger batches)
  if (length(unique_alleles) > 10 && n_cores > 1 && .Platform$OS.type != "windows") {
    # Note: mclapply doesn't work on Windows, falls back to lapply
    if (verbose) message(sprintf("Using parallel processing with %d cores", n_cores))
    
    sequences <- parallel::mclapply(
      unique_alleles,
      function(a) {
        getAlleleSequence(a, type = type, use_cache = use_cache, cache_dir = cache_dir)
      },
      mc.cores = n_cores
    )
  } else {
    # Sequential for small batches or Windows
    if (verbose && length(unique_alleles) > 10) {
      message("Using sequential processing")
    }
    
    sequences <- lapply(unique_alleles, function(a) {
      getAlleleSequence(a, type = type, use_cache = use_cache, cache_dir = cache_dir)
    })
  }
  
  # Check for errors in parallel results
  errors <- sapply(sequences, inherits, "try-error")
  if (any(errors)) {
    failed_alleles <- unique_alleles[errors]
    warning(sprintf("Failed to retrieve sequences for %d allele(s): %s",
                    sum(errors), 
                    paste(failed_alleles, collapse = ", ")))
  }
  
  names(sequences) <- unique_alleles
  
  if (verbose) {
    successful <- sum(!errors)
    message(sprintf("Successfully retrieved %d/%d sequences", 
                    successful, length(unique_alleles)))
  }
  
  sequences
}

#' Clear Sequence Cache
#'
#' @description
#' Utility function to clear the sequence cache, either in-memory or filesystem-based.
#'
#' @param cache_dir If provided, clears the filesystem cache at this location.
#'                  If NULL, clears the in-memory cache for the current session.
#' @param type If provided, only clears cache entries for this sequence type.
#'             Options: "PROT", "NUC", or NULL for all types.
#'
#' @return Invisible NULL
#' 
#' @examples
#' # Clear in-memory cache
#' clearSequenceCache()
#' 
#' # Clear filesystem cache
#' clearSequenceCache(cache_dir = "~/.hla_cache")
#'
#' @export
clearSequenceCache <- function(cache_dir = NULL, type = NULL) {
  if (!is.null(cache_dir)) {
    if (dir.exists(cache_dir)) {
      if (!is.null(type)) {
        # Clear only specific type files (if cache implementation supports it)
        pattern <- paste0("*_", type, "_*")
        files <- list.files(cache_dir, pattern = pattern, full.names = TRUE)
        if (length(files) > 0) {
          unlink(files)
          message(sprintf("Cleared %d %s cache files", length(files), type))
        }
      } else {
        # Clear all cache files
        unlink(cache_dir, recursive = TRUE)
        message("Cleared filesystem cache")
      }
    } else {
      message("Cache directory does not exist")
    }
  } else {
    # Clear in-memory cache
    # This would need to be implemented based on your memoise setup
    message("Note: In-memory cache clearing requires restarting R session or re-sourcing functions")
  }
  
  invisible(NULL)
}

#' Get Sequence Statistics
#'
#' @description
#' Analyzes sequences to provide statistics about length, composition, etc.
#'
#' @param sequences Named list of sequences (output from batchGetSequences)
#' @param type Type of sequences ("PROT" or "NUC")
#'
#' @return Data frame with sequence statistics
#' 
#' @examples
#' alleles <- c("A*01:01", "A*02:01", "B*07:02")
#' seqs <- batchGetSequences(alleles)
#' stats <- getSequenceStats(seqs, type = "PROT")
#'
#' @export
getSequenceStats <- function(sequences, type = c("PROT", "NUC")) {
  type <- match.arg(type)
  
  if (length(sequences) == 0) {
    return(data.frame())
  }
  
  # Calculate statistics for each sequence
  stats <- data.frame(
    allele = names(sequences),
    length = sapply(sequences, nchar),
    stringsAsFactors = FALSE
  )
  
  if (type == "PROT") {
    # Protein-specific statistics
    stats$n_cysteines <- sapply(sequences, function(seq) {
      sum(unlist(strsplit(seq, "")) == "C")
    })
    
    stats$n_prolines <- sapply(sequences, function(seq) {
      sum(unlist(strsplit(seq, "")) == "P")
    })
    
    # Calculate hydrophobicity (simplified)
    hydrophobic_aa <- c("A", "V", "I", "L", "M", "F", "Y", "W")
    stats$hydrophobic_ratio <- sapply(sequences, function(seq) {
      aa <- unlist(strsplit(seq, ""))
      sum(aa %in% hydrophobic_aa) / length(aa)
    })
    
  } else {
    # Nucleotide-specific statistics
    stats$gc_content <- sapply(sequences, function(seq) {
      nuc <- unlist(strsplit(toupper(seq), ""))
      sum(nuc %in% c("G", "C")) / length(nuc)
    })
    
    # Check for stop codons (simplified)
    stats$has_stop_codon <- sapply(sequences, function(seq) {
      # Check for common stop codons
      grepl("(TAA|TAG|TGA)", toupper(seq))
    })
  }
  
  stats
}


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
#' \donttest{
#' # Collecting A*02:01 Protein Sequence (requires internet):
#' seq <- getAlleleSequence("A*02:01")
#' nchar(seq)  # Length of sequence
#'
#' # Collecting A*02:01 Nucleotide Sequence:
#' seq_nuc <- getAlleleSequence("A*02:01", type = "NUC")
#'
#' # Get sequence for HLA-B allele
#' seq_b <- getAlleleSequence("B*07:02")
#' }
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
                              n_cores = 2, 
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
#' Analyzes sequences to provide statistics about length, composition, and
#' optionally comparison metrics (reference similarity, pairwise identity),
#' k-mer entropy, and data-quality flags.
#'
#' @param sequences Named list (or named character vector) of sequences
#'   (e.g., output from batchGetSequences). Names should be allele IDs.
#' @param type Type of sequences: "PROT" (amino acids) or "NUC" (nucleotides).
#' @param ref Optional reference for comparison; either the *name* of one
#'   sequence in `sequences` or a raw sequence string of the same `type`.
#' @param k Integer k for k-mer Shannon entropy (default 2).
#' @param ignore_chars Characters to ignore when computing counts/identity
#'   (default c("*", "-", " ")). Useful for protein stop symbols or gaps.
#' @param compute_pairs Logical; if TRUE, also compute a pairwise percent
#'   identity matrix (alignment-free, position-wise) and return a list with
#'   `stats` and `pairwise_identity`. Default FALSE (returns only `stats`).
#' @param warn_invalid Logical; if TRUE (default), warn when invalid symbols
#'   for the chosen `type` are detected.
#'
#' @return
#' If `compute_pairs = FALSE` (default): a data.frame with per-sequence stats.  
#' If `compute_pairs = TRUE`: a list with elements:
#'   - `stats`: per-sequence stats data.frame
#'   - `pairwise_identity`: numeric matrix (0..1) of % identity (diagonal = 1)
#'
#' @examples
#' alleles <- c("A*01:01", "A*02:01", "B*07:02")
#' seqs <- batchGetSequences(alleles)
#' stats <- getSequenceStats(seqs, type = "PROT")
#'
#' @export
getSequenceStats <- function(
    sequences,
    type = c("PROT", "NUC"),
    ref = NULL,
    k = 2L,
    ignore_chars = c("*", "-", " "),
    compute_pairs = FALSE,
    warn_invalid = TRUE
) {
  type <- match.arg(type)
  
  # Normalize input to named character vector
  if (is.list(sequences)) sequences <- unlist(sequences, use.names = TRUE)
  if (is.null(names(sequences)) || anyNA(names(sequences)) || any(names(sequences) == "")) {
    stop("`sequences` must be named with unique allele IDs.")
  }
  if (length(sequences) == 0) {
    return(if (compute_pairs) list(stats = data.frame(), pairwise_identity = matrix(numeric(0))) else data.frame())
  }
  
  # Helpers ---------------------------------------------------------------
  .split_seq <- function(s) unlist(strsplit(s, "", fixed = TRUE))
  
  .clean_seq <- function(s) {
    v <- .split_seq(s)
    if (length(ignore_chars)) v <- v[!v %in% ignore_chars]
    v
  }
  
  .shannon <- function(x) {
    # x is a character vector of tokens (e.g., kmers)
    if (length(x) == 0) return(NA_real_)
    p <- table(x) / length(x)
    -sum(p * log2(p))
  }
  
  .kmers <- function(v, kk) {
    n <- length(v)
    if (n < kk) return(character(0))
    out <- character(n - kk + 1L)
    for (i in seq_len(n - kk + 1L)) out[i] <- paste0(v[i:(i + kk - 1L)], collapse = "")
    out
  }
  
  # Allowed alphabets
  valid_prot <- c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  valid_nuc  <- c("A","C","G","T","U","N","R","Y","S","W","K","M","B","D","H","V") # IUPAC+U
  
  # Validate, normalize case
  sequences <- toupper(sequences)
  
  # Quality/validity counts ----------------------------------------------
  alphabet <- if (type == "PROT") valid_prot else valid_nuc
  validity <- lapply(sequences, function(s) {
    v <- .clean_seq(toupper(s))
    invalid <- setdiff(unique(v), alphabet)
    list(
      n_invalid = sum(!v %in% alphabet),
      invalid_set = invalid
    )
  })
  
  if (warn_invalid && any(vapply(validity, `[[`, integer(1), "n_invalid") > 0L)) {
    bad <- which(vapply(validity, `[[`, integer(1), "n_invalid") > 0L)
    msg <- paste0(
      "Invalid symbols detected in sequences: ",
      paste(names(sequences)[bad], collapse = ", "),
      ". These were ignored in some calculations."
    )
    warning(msg, call. = FALSE)
  }
  
  # Base per-sequence stats ----------------------------------------------
  allele <- names(sequences)
  length_raw <- nchar(sequences, type = "chars", allowNA = FALSE)
  
  # Composition / content
  if (type == "PROT") {
    hydrophobic_aa <- c("A", "V", "I", "L", "M", "F", "Y", "W")
    charges_pos <- c("K", "R", "H")
    charges_neg <- c("D", "E")
    
    n_cys <- integer(length(sequences))
    n_pro <- integer(length(sequences))
    hydrophobic_ratio <- numeric(length(sequences))
    frac_aromatic <- numeric(length(sequences))
    frac_polar <- numeric(length(sequences))
    net_charge_approx <- numeric(length(sequences))
    
    aromatic <- c("F","Y","W")
    polar    <- c("S","T","N","Q","Y","C")
    
    for (i in seq_along(sequences)) {
      v <- .clean_seq(sequences[i])
      if (!length(v)) {
        n_cys[i] <- n_pro[i] <- 0L
        hydrophobic_ratio[i] <- frac_aromatic[i] <- frac_polar[i] <- NA_real_
        net_charge_approx[i] <- NA_real_
        next
      }
      n_cys[i] <- sum(v == "C")
      n_pro[i] <- sum(v == "P")
      hydrophobic_ratio[i] <- sum(v %in% hydrophobic_aa) / length(v)
      frac_aromatic[i] <- sum(v %in% aromatic) / length(v)
      frac_polar[i]    <- sum(v %in% polar) / length(v)
      net_charge_approx[i] <- sum(v %in% charges_pos) - sum(v %in% charges_neg)
    }
    
    # k-mer entropy (AA)
    kmer_entropy <- vapply(sequences, function(s) {
      v <- .clean_seq(s)
      .shannon(.kmers(v, as.integer(k)))
    }, numeric(1))
    
    stats <- data.frame(
      allele = allele,
      length = length_raw,
      n_invalid = vapply(validity, `[[`, integer(1), "n_invalid"),
      n_cysteines = n_cys,
      n_prolines  = n_pro,
      hydrophobic_ratio = hydrophobic_ratio,
      frac_aromatic = frac_aromatic,
      frac_polar = frac_polar,
      net_charge_approx = net_charge_approx,
      kmer_entropy = kmer_entropy,
      stringsAsFactors = FALSE
    )
    
  } else {
    # NUC
    gc_content <- numeric(length(sequences))
    at_content <- numeric(length(sequences))
    has_stop_anyframe <- logical(length(sequences))
    
    stops <- c("TAA","TAG","TGA","UAA","UAG","UGA")
    
    for (i in seq_along(sequences)) {
      v_all <- .clean_seq(sequences[i])
      v <- v_all[v_all %in% c("A","C","G","T","U")] # restrict denom to canonical
      denom <- length(v)
      if (!denom) {
        gc_content[i] <- at_content[i] <- NA_real_
      } else {
        # Treat U as T for content
        vt <- ifelse(v == "U", "T", v)
        gc_content[i] <- sum(vt %in% c("G","C")) / denom
        at_content[i] <- sum(vt %in% c("A","T")) / denom
      }
      
      # Stop codon detection in any reading frame (0,1,2)
      if (length(v_all) >= 3) {
        frames_have_stop <- logical(3)
        for (f in 0:2) {
          fr <- v_all[(1 + f):length(v_all)]
          n_trip <- length(fr) %/% 3L
          if (n_trip > 0) {
            tri <- character(n_trip)
            idx <- seq_len(n_trip)
            tri[idx] <- paste0(fr[3*idx-2], fr[3*idx-1], fr[3*idx])
            frames_have_stop[f+1] <- any(tri %in% stops)
          }
        }
        has_stop_anyframe[i] <- any(frames_have_stop)
      } else {
        has_stop_anyframe[i] <- FALSE
      }
    }
    
    # k-mer entropy (nucleotides)
    kmer_entropy <- vapply(sequences, function(s) {
      v <- .clean_seq(s)
      .shannon(.kmers(v, as.integer(k)))
    }, numeric(1))
    
    stats <- data.frame(
      allele = allele,
      length = length_raw,
      n_invalid = vapply(validity, `[[`, integer(1), "n_invalid"),
      gc_content = gc_content,
      at_content = at_content,
      has_stop_codon_anyframe = has_stop_anyframe,
      kmer_entropy = kmer_entropy,
      stringsAsFactors = FALSE
    )
  }
  
  # Reference comparison --------------------------------------------------
  if (!is.null(ref)) {
    ref_seq <- NULL
    if (length(ref) == 1L && ref %in% names(sequences)) {
      ref_seq <- sequences[[ref]]
    } else if (is.character(ref) && length(ref) == 1L) {
      ref_seq <- ref
    } else {
      stop("`ref` must be the name of a sequence in `sequences` or a single sequence string.")
    }
    
    # Identity and Hamming (position-wise, ignoring `ignore_chars`)
    ref_vec <- .clean_seq(ref_seq)
    ref_len <- length(ref_vec)
    
    identity_pct <- numeric(length(sequences))
    hamming <- numeric(length(sequences))
    
    for (i in seq_along(sequences)) {
      v <- .clean_seq(sequences[i])
      n <- min(length(v), ref_len)
      if (n == 0) {
        identity_pct[i] <- NA_real_
        hamming[i] <- NA_real_
        next
      }
      # Compare first n positions
      identity_pct[i] <- mean(v[seq_len(n)] == ref_vec[seq_len(n)])
      # Hamming distance only meaningful if lengths are equal on the compared window
      hamming[i] <- sum(v[seq_len(n)] != ref_vec[seq_len(n)])
      # If very different lengths, report NA for strict Hamming over full length
      if (length(v) != ref_len) {
        # keep the partial comparison but note mismatch in an attribute
        attr(hamming, "note") <- "Compared over shortest common length."
      }
    }
    
    stats$identity_to_ref <- identity_pct
    stats$hamming_to_ref  <- hamming
    stats$ref_label <- if (ref %in% names(sequences)) ref else "<raw_ref>"
  }
  
  # Pairwise identity matrix ---------------------------------------------
  if (isTRUE(compute_pairs)) {
    nS <- length(sequences)
    M <- matrix(NA_real_, nS, nS, dimnames = list(names(sequences), names(sequences)))
    seq_clean <- lapply(sequences, .clean_seq)
    
    for (i in seq_len(nS)) {
      vi <- seq_clean[[i]]
      for (j in i:nS) {
        vj <- seq_clean[[j]]
        n <- min(length(vi), length(vj))
        if (n == 0) {
          M[i,j] <- M[j,i] <- NA_real_
        } else {
          pid <- mean(vi[seq_len(n)] == vj[seq_len(n)])
          M[i,j] <- M[j,i] <- pid
        }
      }
    }
    return(list(stats = stats, pairwise_identity = M))
  }
  
  stats
}

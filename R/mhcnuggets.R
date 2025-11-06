#' Predict peptide–MHC binding with mhcnuggets 
#'
#' Calls Python's \code{mhcnuggets.src.predict.predict} inside a basilisk-managed
#' environment to score peptides against a given MHC allele. By default the
#' predictions are written by Python to a temporary CSV and read back in R, which
#' avoids Python stdout capture and is typically faster and more robust.
#'
#' This wrapper exposes most of mhcnuggets' arguments so advanced users can fully
#' control model choice, thresholds, and output format. It also applies a small,
#' session-local patch so older mhcnuggets code that calls
#' \code{keras.optimizers.Adam(lr=...)} works on newer Keras (maps \code{lr}
#' to \code{learning_rate}).
#'
#' @param peptides Character vector of peptide sequences (one per peptide).
#' @param allele MHC allele string. If \code{normalize_allele=TRUE}, common forms
#'   like \code{"A0201"} are normalized to \code{"HLA-A02:01"} for class I and
#'   \code{"HLA-DRB101:01"} for class II.
#' @param mhc_class Either \code{"I"} or \code{"II"} (default \code{"I"}).
#' @param output_path Optional file path for Python to write CSV results. If
#'   \code{NULL} (default), a secure temporary file is used and deleted on exit.
#' @param normalize_allele Logical; normalize \code{allele} to mhcnuggets'
#'   expected format (default \code{TRUE}). Set \code{FALSE} if you already use
#'   exact mhcnuggets allele names.
#' @param model Model architecture string (default \code{"lstm"}).
#' @param mass_spec Logical; use MS-calibrated settings (default \code{FALSE}).
#' @param ic50_threshold Numeric IC50 threshold (nM) for binding calls
#'   (default \code{500}).
#' @param max_ic50 Numeric max IC50 (nM) for capping (default \code{50000}).
#' @param embed_peptides Logical; use embedding (default \code{FALSE}).
#' @param binary_preds Logical; request binary predictions (default \code{FALSE}).
#' @param ba_models Logical; force binding affinity models (default \code{FALSE}).
#' @param rank_output Logical; request rank output (default \code{FALSE}).
#' @param hla_env A \pkg{basilisk} environment object that contains Python +
#'   mhcnuggets (e.g., \code{hlaFerretEnv}).
#'
#' @return A \code{data.frame}. Columns depend on options:
#' \itemize{
#'   \item Default: \code{peptide}, \code{ic50}
#'   \item If \code{binary_preds=TRUE}: \code{peptide}, \code{binary_pred} (plus \code{ic50} if emitted by model)
#'   \item If \code{rank_output=TRUE}: includes \code{rank} (0–1 or percentile)
#' }
#'
#' @examples
#' res <- predictMHCnuggets(
#'   peptides = c("SIINFEKL","LLFGYPVYV"),
#'   allele   = "A*02:01",
#'   mhc_class = "I",
#'   rank_output = TRUE
#' )
#' head(res)
#'
#' @section License and Citation:
#' mhcnuggets is licensed under the GNU General Public License v3.0.
#' If you use MHCnuggets in your work, please cite:
#' Shao, B., et al. (2020). MHCnuggets: A deep learning method for peptide-MHC binding prediction. bioRxiv.
#' GitHub: https://github.com/KarchinLab/mhcnuggets
#'
#' @export
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
#' @importFrom reticulate py_run_string import
#' @importFrom utils read.csv
predictMHCnuggets <- function(peptides,
                              allele,
                              mhc_class        = "I",
                              output_path      = NULL,
                              normalize_allele = TRUE,
                              model            = "lstm",
                              mass_spec        = FALSE,
                              ic50_threshold   = 500,
                              max_ic50         = 50000,
                              embed_peptides   = FALSE,
                              binary_preds     = FALSE,
                              ba_models        = FALSE,
                              rank_output      = FALSE,
                              hla_env          = deepmatchrEnv()) {
  # ---- fast input checks ----
  if (!is.character(peptides)) stop("`peptides` must be character.")
  if (!length(peptides)) return(data.frame(peptide = character(0L)))
  if (!is.character(allele) || length(allele) != 1L) stop("`allele` must be a single string.")
  mhc_class <- toupper(mhc_class)
  if (!(mhc_class %in% c("I","II"))) stop("`mhc_class` must be 'I' or 'II'.")
  
  # ---- optional allele normalization (cheap, common cases) ----
  norm_allele <- function(cls, x) {
    # Uppercase and strip spaces
    a <- toupper(gsub("\\s+", "", x, perl = TRUE))
    # Ensure HLA- prefix
    if (!startsWith(a, "HLA-")) {
      a <- paste0("HLA-", a)
    }
    if (cls == "I") {
      # Convert bare A0201 -> HLA-A02:01
      # Convert HLA-A*02:01 -> HLA-A02:01
      a <- sub("^(HLA-[ABC])\\*?(\\d{2}):(\\d{2})$", "\\1\\2:\\3", a, perl = TRUE)
      a <- sub("^(HLA-[ABC])\\*?(\\d{2})(\\d{2})$", "\\1\\2:\\3", a, perl = TRUE)
    } else {
      # Class II, e.g., DRB1
      # Convert HLA-DRB1*01:01 -> HLA-DRB101:01
      a <- sub("^(HLA-DRB1)\\*?(\\d{2}):(\\d{2})$", "\\1\\2:\\3", a, perl = TRUE)
      a <- sub("^(HLA-DRB1)\\*?(\\d{2})(\\d{2})$", "\\1\\2:\\3", a, perl = TRUE)
    }
    a
  }
  allele_use <- if (isTRUE(normalize_allele)) norm_allele(mhc_class, allele) else allele
  
  # ---- write peptides as a plain text file (one per line) ----
  input_txt <- tempfile("mhcnuggets_peps_", fileext = ".txt")
  # writeLines is efficient and uses minimal copies
  writeLines(peptides, con = input_txt, useBytes = TRUE)
  
  # Output path: user-specified or temp
  tmp_out <- is.null(output_path)
  if (tmp_out) {
    output_path <- tempfile("mhcnuggets_out_", fileext = ".csv")
  }
  # cleanup files on exit
  on.exit({
    try(unlink(input_txt),  silent = TRUE)
    if (tmp_out) try(unlink(output_path), silent = TRUE)
  }, add = TRUE)
  
  # ---- start Python env ----
  proc <- basilisk::basiliskStart(hla_env)
  on.exit(basilisk::basiliskStop(proc), add = TRUE)
  
  # ---- inner Python runner ----
  py_runner <- function(peptides_path, out_path, cls, mhc, model, 
                        model_weights_path, mass_spec, ic50_threshold, 
                        max_ic50, embed_peptides, binary_preds, ba_models, 
                        rank_output) {
    
    # Patch Keras Adam lr->learning_rate to tolerate older mhcnuggets code
    reticulate::py_run_string("
import os, numpy as np
os.environ.setdefault('KERAS_BACKEND', 'tensorflow')

import tensorflow as tf
from tensorflow import keras as tf_keras

# Try standalone keras too
try:
    import keras as standalone_keras
except Exception:
    standalone_keras = None

# Safety: keep float32 policy
try:
    from tensorflow.keras import mixed_precision
    mixed_precision.set_global_policy('float32')
except Exception:
    pass
try:
    tf_keras.backend.set_floatx('float32')
except Exception:
    pass

# ---- Adam(lr=...) -> Adam(learning_rate=...) shim ----
def _adam_shim_factory(Adam):
    def _shim_Adam(*args, **kwargs):
        if 'lr' in kwargs and 'learning_rate' not in kwargs:
            kwargs['learning_rate'] = kwargs.pop('lr')
        return Adam(*args, **kwargs)
    return _shim_Adam

try:
    tf_keras.optimizers.Adam = _adam_shim_factory(tf_keras.optimizers.Adam)
except Exception:
    pass
if standalone_keras is not None:
    try:
        standalone_keras.optimizers.Adam = _adam_shim_factory(standalone_keras.optimizers.Adam)
    except Exception:
        pass
try:
    import mhcnuggets.src.predict as _p
    if hasattr(_p, 'Adam'):
        _p.Adam = tf_keras.optimizers.Adam
except Exception:
    pass

# ---- Helpers to int-cast shapes/indices ----
def _to_int_list(x):
    try:
        if isinstance(x, (list, tuple)):
            return [int(v) if isinstance(v, (float, np.floating)) else v for v in x]
        if hasattr(x, 'numpy'):
            return [int(v) for v in x.numpy().tolist()]
    except Exception:
        pass
    return x

# ---- Slice shims: coerce begin/size/strides to ints ----
if not getattr(tf, '_deepmatchr_slice_shim_v2', False):
    _orig_tf_slice = tf.slice
    def _shim_slice(input_, begin, size, name=None):
        begin = _to_int_list(begin)
        size  = _to_int_list(size)
        return _orig_tf_slice(input_, begin, size, name=name)
    tf.slice = _shim_slice

    # strided_slice is another common path
    _orig_tf_strided_slice = tf.strided_slice
    def _shim_strided_slice(input_, begin, end, strides, begin_mask=0, end_mask=0,
                            ellipsis_mask=0, new_axis_mask=0, shrink_axis_mask=0, name=None):
        begin   = _to_int_list(begin)
        end     = _to_int_list(end)
        strides = _to_int_list(strides)
        return _orig_tf_strided_slice(input_, begin, end, strides,
                                      begin_mask=begin_mask, end_mask=end_mask,
                                      ellipsis_mask=ellipsis_mask, new_axis_mask=new_axis_mask,
                                      shrink_axis_mask=shrink_axis_mask, name=name)
    tf.strided_slice = _shim_strided_slice
    tf._deepmatchr_slice_shim_v2 = True

# ---- Reshape shim (in case shapes pick up floats via math) ----
if not getattr(tf, '_deepmatchr_reshape_shim', False):
    _orig_tf_reshape = tf.reshape
    def _shim_reshape(tensor, shape, name=None):
        shape = _to_int_list(shape)
        return _orig_tf_reshape(tensor, shape, name=name)
    tf.reshape = _shim_reshape
    tf._deepmatchr_reshape_shim = True
")
    
    predict_mod <- reticulate::import("mhcnuggets.src.predict")
    
    # Call predict with explicit named args (no stdout capture)
    predict_mod$predict(
      class_ = cls,
      peptides_path = peptides_path,
      mhc = mhc,
      model = model,
      output = out_path,
      mass_spec = mass_spec,
      ic50_threshold = as.numeric(ic50_threshold),
      max_ic50 = as.numeric(max_ic50),
      embed_peptides = embed_peptides,
      binary_preds = binary_preds,
      ba_models = ba_models,
      rank_output = rank_output
    )
    
    # Return path existence status (cheap boolean)
    file.exists(out_path)
  }
  
  ok <- basilisk::basiliskRun(
    proc, fun = py_runner,
    peptides_path = normalizePath(input_txt, winslash = "/", mustWork = TRUE),
    out_path     = normalizePath(output_path, winslash = "/", mustWork = FALSE),
    cls = mhc_class,
    mhc = allele_use,
    model = model,
    mass_spec = mass_spec,
    ic50_threshold = ic50_threshold,
    max_ic50 = max_ic50,
    embed_peptides = embed_peptides,
    binary_preds = binary_preds,
    ba_models = ba_models,
    rank_output = rank_output
  )
  
  if (!isTRUE(ok) || !file.exists(output_path)) {
    stop("mhcnuggets did not produce an output file; check allele '", allele_use,
         "', model weights, and environment configuration.")
  }
  
  # ---- parse CSV with minimal overhead ----
  df <- utils::read.csv(output_path, stringsAsFactors = FALSE)
  
  # Harmonize common variants
  if (!("peptide" %in% names(df))) {
    stop("Unexpected mhcnuggets output: missing 'peptide' column.")
  }
  # If ic50 named differently (rare), try to coerce
  if (!("ic50" %in% names(df))) {
    if ("affinity" %in% names(df)) names(df)[names(df)=="affinity"] <- "ic50"
  }
  # Ensure numeric where present
  if ("ic50" %in% names(df)) df$ic50 <- suppressWarnings(as.numeric(df$ic50))
  if ("rank" %in% names(df)) df$rank <- suppressWarnings(as.numeric(df$rank))
  
  df
}

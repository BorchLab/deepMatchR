# Linux/x86_64 (Bioconductor builders)
.deepmatchrEnv_linux <- basilisk::BasiliskEnvironment(
  envname = "deepmatchrEnv_v2",   # bump name to force a fresh env on CI
  pkgname = "deepMatchR",
  packages = c("python=3.10"),
  pip = c(
    "mhcnuggets==2.4.1",
    "tensorflow==2.15.1",    # Linux CPU build
    "keras==2.15.0",         # optional; OK to omit since TF 2.15 vendors Keras
    "numpy==1.26.4",         # instead of "numpy<2.0"
    "protobuf==3.20.3"       # instead of "protobuf<4"
  )
)

# macOS Apple Silicon (developers on M-series Macs)
.deepmatchrEnv_macos <- basilisk::BasiliskEnvironment(
  envname = "deepmatchrEnv_v2",   # keep the same bumped name
  pkgname = "deepMatchR",
  packages = c("python=3.10"),
  pip = c(
    "mhcnuggets==2.4.1",
    "tensorflow-macos==2.15.0",
    "keras==2.15.0",
    "numpy==1.26.4",
    "protobuf==3.20.3"
  )
)

#' Return a Basilisk environment tailored to the current platform
#' @param platform "auto", "linux", or "macos"
#' @export
deepmatchrEnv <- function(platform = c("auto","linux","macos")) {
  platform <- match.arg(platform)
  if (platform == "auto") {
    os <- tolower(Sys.info()[["sysname"]] %||% .Platform$OS.type)
    if (grepl("darwin|mac", os)) return(.deepmatchrEnv_macos)
    return(.deepmatchrEnv_linux)
  }
  if (platform == "macos") return(.deepmatchrEnv_macos)
  .deepmatchrEnv_linux
}

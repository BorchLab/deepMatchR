# Linux/x86_64 (Bioconductor builders)
.deepmatchrEnv_linux <- basilisk::BasiliskEnvironment(
  envname = "deepmatchrEnv",
  pkgname = "deepMatchR",
  packages = c("python=3.10"),
  pip = c(
    "mhcnuggets==2.4.1",
    "tensorflow==2.15.1",  
    "keras==2.15.0",
    "numpy<2.0",          
    "protobuf<4"           
  )
)


# macOS Apple Silicon (for developers/users on M-series Macs)
.deepmatchrEnv_macos <- basilisk::BasiliskEnvironment(
  envname = "deepmatchrEnv",
  pkgname = "deepMatchR",
  packages = c("python=3.10"),
  pip = c(
    "mhcnuggets==2.4.1",
    "tensorflow-macos==2.15.0",
    "keras==2.15.0",
    "numpy<2.0",
    "protobuf<4"
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
    return(.deepmatchrEnv_linux) # default to Linux
  }
  if (platform == "macos") return(.deepmatchrEnv_macos)
  .deepmatchrEnv_linux
}
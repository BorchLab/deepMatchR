# Linux/x86_64 (Bioconductor builders)
.deepmatchrEnv_linux <- basilisk::BasiliskEnvironment(
  envname = "deepmatchrEnv_v2",   # bump name to force a fresh env on CI
  pkgname = "deepMatchR",
  packages = c("python=3.10"),
  pip = c(
    "mhcnuggets==2.4.1",
    "tensorflow==2.19.1"   # Linux CPU build
  )
)

# macOS Apple Silicon (developers on M-series Macs)
.deepmatchrEnv_macos <- basilisk::BasiliskEnvironment(
  envname = "deepmatchrEnv_v2",   # keep the same bumped name
  pkgname = "deepMatchR",
  packages = c("python=3.10"),
  pip = c(
    "mhcnuggets==2.4.1",
    "tensorflow-macos==2.16.1"
  )
)

#' Return a Basilisk environment tailored to the current platform
#'
#' @description
#' Returns the appropriate Basilisk environment for MHCnuggets based on
#' the current operating system. This environment is used internally by
#' \code{\link{predictMHCnuggets}} for Python integration.
#'
#' @param platform Character. One of "auto" (default), "linux", or "macos".
#'   When "auto", the platform is detected automatically.
#'
#' @return A \code{BasiliskEnvironment} object configured for the current platform.
#'
#' @examples
#' # Get the environment for the current platform
#' env <- deepmatchrEnv()
#' print(class(env))
#'
#' # Explicitly request Linux environment
#' env_linux <- deepmatchrEnv(platform = "linux")
#'
#' @seealso \code{\link{predictMHCnuggets}}
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

# deploy.R - Build and deploy pkgdown site to borch.dev
#
# Usage:
#   source("pkgdown/deploy.R")
#   # or from command line:
#   Rscript pkgdown/deploy.R
#
# This script:
# 1. Builds the pkgdown site locally
# 2. Copies it to the borch.dev website folder
# 3. Optionally commits changes to the website repo

# Configuration
WEBSITE_PATH <- "/Users/nick/Documents/GitHub/borcherding/static/uploads/deepMatchR"
PACKAGE_ROOT <- here::here()

# Build the site
message("Building pkgdown site...")
pkgdown::build_site(pkg = PACKAGE_ROOT, preview = FALSE)

# Get the docs folder
docs_path <- file.path(PACKAGE_ROOT, "docs")

if (!dir.exists(docs_path)) {
 stop("docs/ folder not found. Did pkgdown::build_site() fail?")
}

# Create destination if it doesn't exist
if (!dir.exists(WEBSITE_PATH)) {
  message("Creating destination folder: ", WEBSITE_PATH)
  dir.create(WEBSITE_PATH, recursive = TRUE)
}

# Sync files (remove old, copy new)
message("Syncing to website folder...")
message("  From: ", docs_path)
message("  To:   ", WEBSITE_PATH)

# Remove old files (but preserve .DS_Store and .git if any)
old_files <- list.files(WEBSITE_PATH, full.names = TRUE, all.files = FALSE)
if (length(old_files) > 0) {
  unlink(old_files, recursive = TRUE)
}

# Copy new files
file.copy(
  from = list.files(docs_path, full.names = TRUE, all.files = TRUE, include.dirs = TRUE),
  to = WEBSITE_PATH,
  recursive = TRUE,
  overwrite = TRUE
)
message("Done! Files copied to: ", WEBSITE_PATH)

# Optional: Show git status of website repo
website_repo <- dirname(dirname(dirname(WEBSITE_PATH)))
message("\nWebsite repo status:")
message("  Path: ", website_repo)
tryCatch({
  system2("git", args = c("-C", website_repo, "status", "--short", "static/uploads/deepMatchR"),
          stdout = TRUE, stderr = TRUE) |> cat(sep = "\n")
}, error = function(e) {
  message("  (Could not check git status)")
})

message("\nTo commit changes to website repo:")
message('  cd "', website_repo, '"')
message('  git add static/uploads/deepMatchR')
message('  git commit -m "Update deepMatchR pkgdown site"')
message('  git push')

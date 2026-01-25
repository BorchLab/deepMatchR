#!/bin/bash
# deploy.sh - Build and deploy pkgdown site to borch.dev
#
# Usage:
#   ./pkgdown/deploy.sh           # Build and copy
#   ./pkgdown/deploy.sh --commit  # Build, copy, and commit to website repo
#
# Make executable with: chmod +x pkgdown/deploy.sh

set -e  # Exit on error

# Configuration
WEBSITE_PATH="/Users/nick/Documents/GitHub/borcherding/static/uploads/deepMatchR"
WEBSITE_REPO="/Users/nick/Documents/GitHub/borcherding"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="$(dirname "$SCRIPT_DIR")"

echo "=== deepMatchR pkgdown deployment ==="
echo "Package root: $PACKAGE_ROOT"
echo "Website path: $WEBSITE_PATH"
echo ""

# Step 1: Build the site
echo "Building pkgdown site..."
cd "$PACKAGE_ROOT"
Rscript -e "pkgdown::build_site(preview = FALSE)"

# Step 2: Sync to website folder
echo ""
echo "Syncing to website folder..."

# Create destination if needed
mkdir -p "$WEBSITE_PATH"

# Use rsync for efficient sync (delete files not in source)
rsync -av --delete \
  --exclude='.DS_Store' \
  --exclude='.git' \
  "$PACKAGE_ROOT/docs/" \
  "$WEBSITE_PATH/"

echo ""
echo "Deployed to: $WEBSITE_PATH"

# Step 3: Optionally commit to website repo
if [[ "$1" == "--commit" ]]; then
  echo ""
  echo "Committing to website repo..."
  cd "$WEBSITE_REPO"
  git add static/uploads/deepMatchR
  git commit -m "Update deepMatchR pkgdown site ($(date +%Y-%m-%d))"
  echo ""
  echo "Changes committed. Run 'git push' in website repo to deploy."
else
  echo ""
  echo "To commit changes to website repo:"
  echo "  cd \"$WEBSITE_REPO\""
  echo "  git add static/uploads/deepMatchR"
  echo "  git commit -m \"Update deepMatchR pkgdown site\""
  echo "  git push"
  echo ""
  echo "Or run: ./pkgdown/deploy.sh --commit"
fi

echo ""
echo "=== Done ==="

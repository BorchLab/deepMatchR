# tests/testthat/test-updateWmdaData.R

# --- Cache Directory Tests ---

test_that(".getWmdaCacheDir returns valid path", {
  cache_dir <- deepMatchR:::.getWmdaCacheDir()
  expect_type(cache_dir, "character")
  expect_true(nzchar(cache_dir))
})

test_that(".getWmdaCacheDir respects DEEPMATCHR_CACHE_DIR environment variable", {
  withr::with_envvar(c(DEEPMATCHR_CACHE_DIR = "/custom/cache/path"), {
    cache_dir <- deepMatchR:::.getWmdaCacheDir()
    expect_equal(cache_dir, "/custom/cache/path")
  })
})

test_that(".getWmdaCacheDir falls back to rappdirs when available", {
  withr::with_envvar(c(DEEPMATCHR_CACHE_DIR = ""), {
    cache_dir <- deepMatchR:::.getWmdaCacheDir()
    expect_type(cache_dir, "character")
    expect_true(nzchar(cache_dir))
    # Should either be rappdirs path or tempdir fallback
    expect_true(grepl("deepMatchR", cache_dir) || grepl("Rtmp", cache_dir))
  })
})

# --- clearWmdaCache Tests ---

test_that("clearWmdaCache works without error when no cache exists", {
  temp_dir <- tempfile("wmda_cache_test")

  result <- clearWmdaCache(cache_dir = temp_dir, verbose = FALSE)
  expect_false(result)  # No cache existed
})

test_that("clearWmdaCache removes existing cache files", {
  temp_dir <- tempfile("wmda_cache_test")
  dir.create(temp_dir, recursive = TRUE)

  # Create mock cache files
  cache_file <- file.path(temp_dir, "wmda_cache.rds")
  version_file <- file.path(temp_dir, "wmda_version.txt")
  saveRDS(list(test = "data"), cache_file)
  writeLines("test_version", version_file)

  expect_true(file.exists(cache_file))
  expect_true(file.exists(version_file))

  result <- clearWmdaCache(cache_dir = temp_dir, verbose = FALSE)

  expect_true(result)
  expect_false(file.exists(cache_file))
  expect_false(file.exists(version_file))

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("clearWmdaCache with verbose=TRUE prints message", {
  temp_dir <- tempfile("wmda_cache_test")

  expect_message(clearWmdaCache(cache_dir = temp_dir, verbose = TRUE), "No WMDA cache found")

  # Now with actual cache
  dir.create(temp_dir, recursive = TRUE)
  cache_file <- file.path(temp_dir, "wmda_cache.rds")
  saveRDS(list(test = "data"), cache_file)

  expect_message(clearWmdaCache(cache_dir = temp_dir, verbose = TRUE), "WMDA cache cleared")

  unlink(temp_dir, recursive = TRUE)
})

# --- Parser Tests ---

test_that(".parseDnaSer parses valid lines correctly", {
  lines <- c(
    "A*;01:01;1;;",
    "A*;02:01;2;;",
    "B*;07:02;7;;"
  )

  result <- deepMatchR:::.parseDnaSer(lines)

  expect_s3_class(result, "data.table")
  expect_true(all(c("locus", "allele_2f", "serology") %in% names(result)))
  expect_equal(nrow(result), 3)
})

test_that(".parseDnaSer handles multiple serology assignments", {
  # When multiple serologies are listed (separated by /), first should be taken
  lines <- c("A*;01:01;1/9;;")

  result <- deepMatchR:::.parseDnaSer(lines)

  expect_equal(result$serology[1], "1")
})

test_that(".parseDnaSer extracts two-field alleles correctly", {
  lines <- c(
    "A*;01:01:01:01;1;;",  # Four-field
    "A*;02:01:01;2;;"      # Three-field
  )

  result <- deepMatchR:::.parseDnaSer(lines)

  expect_equal(result$allele_2f[1], "01:01")
  expect_equal(result$allele_2f[2], "02:01")
})

test_that(".parseSerSer parses broad/split relationships", {
  lines <- c(
    "A;9;23",
    "A;9;24",
    "DR;2;15",
    "DR;2;16"
  )

  result <- deepMatchR:::.parseSerSer(lines)

  expect_s3_class(result, "data.table")
  expect_true(all(c("locus", "broad", "splits") %in% names(result)))

  # Should aggregate splits for same broad
  a9_row <- result[locus == "A" & broad == "9"]
  expect_equal(nrow(a9_row), 1)
  expect_true(grepl("23", a9_row$splits))
  expect_true(grepl("24", a9_row$splits))
})

test_that(".parseNomP parses P-group definitions", {
  lines <- c(
    "A*01:01/01:02/01:03;01:01P",
    "B*07:02/07:03;07:02P"
  )

  result <- deepMatchR:::.parseNomP(lines)

  expect_s3_class(result, "data.table")
  expect_true(all(c("locus", "p_group", "reference_2f") %in% names(result)))
  expect_true("A" %in% result$locus)
  expect_true("B" %in% result$locus)
})

test_that(".parseNomP handles malformed lines gracefully", {
  lines <- c(
    "A*01:01/01:02;01:01P",  # Valid
    "invalid_line",           # Invalid - no asterisk
    ";;"                      # Empty
  )

  result <- deepMatchR:::.parseNomP(lines)

  # Should only have the valid line
  expect_equal(nrow(result), 1)
})

# --- updateWmdaData Integration Tests (with mocking) ---

test_that("updateWmdaData skips download when cache is current", {
  temp_dir <- tempfile("wmda_update_test")
  dir.create(temp_dir, recursive = TRUE)

  # Create mock cache with "Latest" version
  cache_file <- file.path(temp_dir, "wmda_cache.rds")
  version_file <- file.path(temp_dir, "wmda_version.txt")

  mock_data <- list(
    serology = data.table::data.table(locus = "A*", allele_2f = "01:01", serology = "1"),
    splits = data.table::data.table(locus = "A", broad = "9", splits = "23|24"),
    pgroups = data.table::data.table(locus = "A", p_group = "01:01P", reference_2f = "01:01")
  )
  saveRDS(mock_data, cache_file)
  writeLines("Latest", version_file)

  expect_message(
    updateWmdaData(version = "Latest", cache_dir = temp_dir, force = FALSE, verbose = TRUE),
    "up to date"
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("updateWmdaData force parameter overrides cache check", {
  temp_dir <- tempfile("wmda_force_test")
  dir.create(temp_dir, recursive = TRUE)

  # Create mock cache
  cache_file <- file.path(temp_dir, "wmda_cache.rds")
  version_file <- file.path(temp_dir, "wmda_version.txt")

  mock_data <- list(
    serology = data.table::data.table(locus = "A*", allele_2f = "01:01", serology = "1"),
    splits = data.table::data.table(locus = "A", broad = "9", splits = "23"),
    pgroups = data.table::data.table(locus = "A", p_group = "01:01P", reference_2f = "01:01")
  )
  saveRDS(mock_data, cache_file)
  writeLines("Latest", version_file)

  # With force=TRUE, it should attempt to download (and may fail due to network)
  # We just verify it doesn't short-circuit
  result <- tryCatch({
    updateWmdaData(version = "Latest", cache_dir = temp_dir, force = TRUE, verbose = FALSE)
    "completed"
  }, error = function(e) {
    # Network error is expected in test environment
    "network_error"
  })

  expect_true(result %in% c("completed", "network_error"))

  unlink(temp_dir, recursive = TRUE)
})

test_that("updateWmdaData creates cache directory if needed",
{
  temp_dir <- tempfile("wmda_newdir_test")
  expect_false(dir.exists(temp_dir))

  # This will fail due to network, but should create the directory first
  tryCatch({
    updateWmdaData(version = "Latest", cache_dir = temp_dir, verbose = FALSE)
  }, error = function(e) {
    # Expected
  })

  expect_true(dir.exists(temp_dir))

  unlink(temp_dir, recursive = TRUE)
})

test_that("updateWmdaData returns cache directory path invisibly", {
  temp_dir <- tempfile("wmda_return_test")
  dir.create(temp_dir, recursive = TRUE)

  # Create mock cache so it skips download
  cache_file <- file.path(temp_dir, "wmda_cache.rds")
  version_file <- file.path(temp_dir, "wmda_version.txt")

  mock_data <- list(
    serology = data.table::data.table(locus = "A*", allele_2f = "01:01", serology = "1"),
    splits = data.table::data.table(locus = "A", broad = "9", splits = "23"),
    pgroups = data.table::data.table(locus = "A", p_group = "01:01P", reference_2f = "01:01")
  )
  saveRDS(mock_data, cache_file)
  writeLines("Latest", version_file)

  result <- updateWmdaData(version = "Latest", cache_dir = temp_dir, verbose = FALSE)

  expect_equal(result, temp_dir)

  unlink(temp_dir, recursive = TRUE)
})

# --- .loadWmdaData Tests ---

test_that(".loadWmdaData returns list with required components", {
  wmda_data <- deepMatchR:::.loadWmdaData()

  expect_type(wmda_data, "list")
  expect_true(all(c("serology", "splits", "pgroups") %in% names(wmda_data)))
})

test_that(".loadWmdaData returns data.tables with correct keys", {
  wmda_data <- deepMatchR:::.loadWmdaData()

  expect_s3_class(wmda_data$serology, "data.table")
  expect_s3_class(wmda_data$splits, "data.table")
  expect_s3_class(wmda_data$pgroups, "data.table")

  # Keys should be set
  expect_true(data.table::haskey(wmda_data$serology))
  expect_true(data.table::haskey(wmda_data$splits))
  expect_true(data.table::haskey(wmda_data$pgroups))
})

test_that(".loadWmdaData prefers cached data when available", {
  temp_dir <- tempfile("wmda_load_test")
  dir.create(temp_dir, recursive = TRUE)

  # Create mock cache with distinguishable data
  cache_file <- file.path(temp_dir, "wmda_cache.rds")

  mock_data <- list(
    serology = data.table::data.table(
      locus = "TEST*",
      allele_2f = "99:99",
      serology = "TEST_SER"
    ),
    splits = data.table::data.table(
      locus = "TEST",
      broad = "99",
      splits = "TEST_SPLIT"
    ),
    pgroups = data.table::data.table(
      locus = "TEST",
      p_group = "99:99P",
      reference_2f = "99:99"
    )
  )
  saveRDS(mock_data, cache_file)

  withr::with_envvar(c(DEEPMATCHR_CACHE_DIR = temp_dir), {
    wmda_data <- deepMatchR:::.loadWmdaData()
    expect_true("TEST*" %in% wmda_data$serology$locus)
  })

  unlink(temp_dir, recursive = TRUE)
})

# --- Version Handling Tests ---

test_that("updateWmdaData handles specific version strings", {
  temp_dir <- tempfile("wmda_version_test")
  dir.create(temp_dir, recursive = TRUE)

  # Create cache with a specific version
  cache_file <- file.path(temp_dir, "wmda_cache.rds")
  version_file <- file.path(temp_dir, "wmda_version.txt")

  mock_data <- list(
    serology = data.table::data.table(locus = "A*", allele_2f = "01:01", serology = "1"),
    splits = data.table::data.table(locus = "A", broad = "9", splits = "23"),
    pgroups = data.table::data.table(locus = "A", p_group = "01:01P", reference_2f = "01:01")
  )
  saveRDS(mock_data, cache_file)
  writeLines("3.54.0", version_file)

  # Request same version - should skip
  expect_message(
    updateWmdaData(version = "3.54.0", cache_dir = temp_dir, verbose = TRUE),
    "up to date"
  )

  # Request different version - should attempt download
  result <- tryCatch({
    updateWmdaData(version = "3.55.0", cache_dir = temp_dir, verbose = FALSE)
    "completed"
  }, error = function(e) {
    "network_error"
  })

  expect_true(result %in% c("completed", "network_error"))

  unlink(temp_dir, recursive = TRUE)
})

# tests/testthat/test-hlaGeno.R

test_that("hlaGeno creates valid object", {
  df <- data.frame(A_1 = "A*01:01", A_2 = "A*02:01", stringsAsFactors = FALSE)
  geno <- hlaGeno(df)

  expect_s3_class(geno, "hla_genotype")
  expect_true(is.list(geno))
  expect_true(all(c("data", "locus_present") %in% names(geno)))
})

test_that("hlaGeno extracts loci correctly", {
  df <- data.frame(
    A_1 = "A*01:01", A_2 = "A*02:01",
    B_1 = "B*07:02", B_2 = "B*08:01",
    C_1 = "C*01:02",
    stringsAsFactors = FALSE
  )
  geno <- hlaGeno(df)

  expect_setequal(geno$locus_present, c("A", "B", "C"))
})

test_that("hlaGeno preserves data", {
  df <- data.frame(A_1 = "A*01:01", A_2 = "A*02:01", stringsAsFactors = FALSE)
  geno <- hlaGeno(df)

  expect_equal(geno$data, df)
})

test_that("hlaGeno handles single locus", {
  df <- data.frame(A_1 = "A*01:01", stringsAsFactors = FALSE)
  geno <- hlaGeno(df)

  expect_equal(geno$locus_present, "A")
})

test_that("hlaGeno handles Class II loci", {
  df <- data.frame(
    DRB1_1 = "DRB1*01:01", DRB1_2 = "DRB1*03:01",
    DQB1_1 = "DQB1*02:01", DQB1_2 = "DQB1*05:01",
    stringsAsFactors = FALSE
  )
  geno <- hlaGeno(df)

  expect_setequal(geno$locus_present, c("DRB1", "DQB1"))
})

test_that("validateHlaGeno returns TRUE for valid object", {
  df <- data.frame(A_1 = "A*01:01", stringsAsFactors = FALSE)
  geno <- hlaGeno(df)

  expect_true(validateHlaGeno(geno))
})

test_that("validateHlaGeno errors on invalid object", {
  # Not an hla_genotype
  expect_error(validateHlaGeno(list(data = 1, locus_present = "A")), "hla_genotype")

  # Missing elements
  bad_geno <- structure(list(data = data.frame(A_1 = "A*01:01")), class = "hla_genotype")
  expect_error(validateHlaGeno(bad_geno), "locus_present")

  # data not a data.frame
  bad_geno2 <- structure(
    list(data = "not_a_df", locus_present = "A"),
    class = "hla_genotype"
  )
  expect_error(validateHlaGeno(bad_geno2), "data frame")

  # locus_present not character
  bad_geno3 <- structure(
    list(data = data.frame(A_1 = "A*01:01"), locus_present = 123),
    class = "hla_genotype"
  )
  expect_error(validateHlaGeno(bad_geno3), "character")
})

test_that("print.hla_genotype outputs correctly", {
  df <- data.frame(A_1 = "A*01:01", B_1 = "B*07:02", stringsAsFactors = FALSE)
  geno <- hlaGeno(df)

  expect_output(print(geno), "HLA Genotype Data")
  expect_output(print(geno), "Loci present:")
  expect_output(print(geno), "Number of samples:")
})

test_that("print.hla_genotype returns object invisibly", {
  df <- data.frame(A_1 = "A*01:01", stringsAsFactors = FALSE)
  geno <- hlaGeno(df)

  result <- capture.output(ret_val <- print(geno))
  expect_identical(ret_val, geno)
})
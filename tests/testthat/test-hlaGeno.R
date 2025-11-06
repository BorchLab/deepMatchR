#tests/testthat/test-hlaGeno.R

test_that("hlaGeno builds a well-formed hla_genotype object", {
  df <- data.frame(
    A_1 = "A*01:01", A_2 = "A*02:01",
    B_1 = "B*07:02", B_2 = "B*08:01",
    stringsAsFactors = FALSE
  )
  
  g <- hlaGeno(df)
  expect_s3_class(g, "hla_genotype")
  expect_true(is.data.frame(g$data))
  expect_setequal(g$locus_present, c("A","B"))
  expect_true(validateHlaGeno(g))
})

test_that("validateHlaGeno rejects malformed objects", {
  expect_error(validateHlaGeno(list()), "hla_genotype")
  bad1 <- structure(list(data = 1, locus_present = c("A")), class = "hla_genotype")
  expect_error(validateHlaGeno(bad1), "data.*data frame")
  bad2 <- structure(list(data = data.frame(), locus_present = 1L), class = "hla_genotype")
  expect_error(validateHlaGeno(bad2), "locus_present.*character")
})

test_that("print.hla_genotype produces expected header and returns invisibly", {
  df <- data.frame(A_1 = "A*01:01", A_2 = "A*02:01")
  g  <- hlaGeno(df)
  
  out <- paste(capture.output(print(g)), collapse = "\n")
  expect_match(out, "HLA Genotype Data")
  expect_match(out, "Loci present:")
  expect_match(out, "Number of samples:")
})

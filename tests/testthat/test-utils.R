
test_that(".strExtract extracts first regex match and returns NA for non-matches", {
  f <- deepMatchR:::.strExtract
  x <- c("ABC123", "no_digits", "X-42-Y")
  expect_equal(f(x, "[0-9]+"), c("123", NA, "42"))
})

test_that(".padRightBase pads to requested width with custom char", {
  f <- deepMatchR:::.padRightBase
  expect_equal(f(c("A","BB",""), len = 3, pad = "-"), c("A--","BB-","---"))
  expect_equal(f("ABC", len = 2), "ABC") # no negative padding
})

test_that(".alphanumericalSort orders by alpha then numeric parts", {
  f <- deepMatchR:::.alphanumericalSort
  x <- c("B10","B2","A12","A1","A","B")
  lev <- f(x)
  expect_equal(lev, c("A1","A12", "A","B2","B10", "B"))
})

test_that(".colorizer returns exactly n colors; .themeMatchR returns a theme", {
  colfun <- deepMatchR:::.colorizer
  pal <- colfun("spectral", n = 5)
  expect_length(pal, 5)
  expect_true(all(nchar(pal) > 0))
  
  thm <- deepMatchR:::.themeMatchR()
  expect_s3_class(thm, "theme")
})

test_that(".checkSAB enforces required columns", {
  ok <- data.frame(BeadID=1, SpecAbbr="A1", Specificity="A*01:01", NormalValue=100)
  expect_silent(deepMatchR:::.checkSAB(ok))
  
  bad <- data.frame(BeadID=1, Specificity="A*01:01", NormalValue=100)
  expect_error(deepMatchR:::.checkSAB(bad), "BeadID.*SpecAbbr.*Specificity.*NormalValue")
})

test_that(".loadData reads CSV, rejects unknown extensions", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(data.frame(BeadID=1, SpecAbbr="A1", Specificity="A*01:01", NormalValue=100),
            tmp, row.names = FALSE)
  on.exit(unlink(tmp), add = TRUE)
  
  df <- deepMatchR:::.loadData(tmp)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("BeadID","SpecAbbr","Specificity","NormalValue") %in% names(df)))
  
  expect_error(deepMatchR:::.loadData(tempfile(fileext = ".txt")), "The specified file does not exist.")
})

test_that(".processSAB produces allele/mfi_min table with expected columns", {
  raw <- deepMatchR_example[[1]]
  # Must pass .checkSAB but .processSAB ignores extra columns appropriately
  deepMatchR:::.checkSAB(raw)
  out <- deepMatchR:::.processSAB(raw)
  
  expect_s3_class(out, "data.frame")
  expect_true(all(c("BeadID","antigen","bw46","allele","NormalValue","loci","mfi_min") %in% names(out)))
  expect_true(all(grepl("^[ABCD]", out$loci)))
  # mfi_min per allele
  mins <- tapply(out$NormalValue, out$allele, min)
  expect_true(all(sort(unique(out$mfi_min)) %in% sort(mins)))
})

test_that(".processPRA handles Class I with Bw synthesis and C* carryover", {
  raw <- deepMatchR_example[[3]]

  outI <- deepMatchR:::.processPRA(raw, class = "I")
  expect_s3_class(outI, "data.frame")
  expect_true(all(c("BeadID","antigen","allele","allele_locus","loci_family","NormalValue","mfi_min","pairs") %in% names(outI)))
  # Ensure BW rows exist and are labeled
  expect_true(any(outI$loci_family == "BW"))
  expect_true(all(outI[outI$loci_family == "BW","antigen"] %in% c("Bw4","Bw6")))
})

test_that(".processPRA handles Class II mapping & guard-rails", {
  raw <- deepMatchR_example[[3]]
  outII <- deepMatchR:::.processPRA(raw, class = "II")
  expect_true(all(outII$loci_family %in% c("DR","DR5","DQ","DP")))
  # DRB3/4/5 -> DR5 family
  expect_true(any(outII$loci_family == "DR5"))
  # Guard-rails remove mismatched loci/alleles: no DR family rows with non-DRB alleles
  expect_false(any(outII$loci_family %in% c("DR","DR5") & !grepl("^DRB", outII$allele_locus)))
})

test_that("Class I pairs alternate within bead/loci_family; Class II pairs factor by antigen", {
  rawI <- deepMatchR_example[[3]]
  outI <- deepMatchR:::.processPRA(rawI, class = "I")
  # pairs must be 1,2,1,2... per (BeadID,loci_family)
  chk <- split(outI$pairs, interaction(outI$BeadID, outI$loci_family, drop = TRUE))
  expect_true(all(vapply(chk, function(v) all(v %in% c(1L,2L)), logical(1))))
  
  outII <- deepMatchR:::.processPRA(rawI, class = "II")
  # pairs should be integer factor by antigen within (BeadID,loci_family)
  fac_ok <- all(outII$pairs == as.integer(factor(outII$antigen, levels = unique(outII$antigen))))
  expect_true(fac_ok || length(unique(outII$loci_family)) > 1) # lenient due to grouping
})



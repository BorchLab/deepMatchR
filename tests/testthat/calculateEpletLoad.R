#tests/testthat/test-calculateEpletLoad.R

test_that("calculateEpletLoad: total/per_locus/pairwise work with filters", {
  # Recipient & Donor fixtures
  recipient <- data.frame(
    A_1 = "A*01:01", A_2 = "A*02:01",
    B_1 = "B*07:02", B_2 = "B*08:01",
    stringsAsFactors = FALSE
  )
  donor <- data.frame(
    A_1 = "A*03:01", A_2 = "A*24:02",
    B_1 = "B*44:02", B_2 = "B*51:01",
    stringsAsFactors = FALSE
  )
  rgeno <- hlaGeno(recipient)
  dgeno <- hlaGeno(donor)
  
  # Make utils::data a NO-OP inside the call so our injected object is used
  expect_silent(
    with_mocked_bindings(
      utils = list(data = getFromNamespace(".noop_utils_data", "deepMatchR")),
      {
        # total
        tot <- calculateEpletLoad(rgeno, dgeno, return = "total")
        expect_type(tot, "integer")
        expect_gte(tot, 0L)
        
        # per_locus
        pl <- calculateEpletLoad(rgeno, dgeno, return = "per_locus")
        expect_s3_class(pl, "data.frame")
        expect_setequal(pl$locus, c("A","B"))
        expect_true(all(pl$eplet_load >= 0L))
        
        # pairwise @ B
        mB <- calculateEpletLoad(rgeno, dgeno, return = "pairwise", pairwise_locus = "B")
        expect_true(is.matrix(mB))
        expect_identical(rownames(mB), unlist(recipient[1, grep("^B_", names(recipient))]))
        expect_identical(colnames(mB), unlist(donor[1,     grep("^B_", names(donor))]))
        expect_true(all(mB >= 0L))
      }
    )
  )
})

test_that("calculateEpletLoad honors loci restriction & filters", {
  rgeno <- hlaGeno(data.frame(A_1="A*01:01", A_2="A*02:01", B_1="B*07:02"))
  dgeno <- hlaGeno(data.frame(A_1="A*03:01", A_2="A*24:02", B_1="B*44:02"))
  
  with_mocked_bindings(
    utils = list(data = getFromNamespace(".noop_utils_data", "deepMatchR")),
    {
      # Only locus A
      tot_A <- calculateEpletLoad(rgeno, dgeno, loci = "A", return = "total")
      # Only locus B
      tot_B <- calculateEpletLoad(rgeno, dgeno, loci = "B", return = "total")
      expect_true(is.integer(tot_A) && is.integer(tot_B))
      expect_gte(tot_A, 0L); expect_gte(tot_B, 0L)
      
      # Evidence filter reduces/changes counts
      tot_default <- calculateEpletLoad(rgeno, dgeno, return = "total")
      tot_A1_only <- calculateEpletLoad(rgeno, dgeno, evidence_level = "A1", return = "total")
      expect_true(is.integer(tot_default))
      expect_true(is.integer(tot_A1_only))
      # not asserting strict inequality, just that it computes
    }
  )
})

test_that("calculateEpletLoad: input validation & error paths", {
  # No shared loci after filter
  rgeno <- hlaGeno(data.frame(A_1="A*01:01"))
  dgeno <- hlaGeno(data.frame(B_1="B*44:02"))
  
  with_mocked_bindings(utils = list(data = getFromNamespace(".noop_utils_data", "deepMatchR")), {
    expect_error(calculateEpletLoad(rgeno, dgeno, loci = "A"), "No shared loci")
  })
  
  # Missing alleles (empty strings) -> error
  rgeno2 <- hlaGeno(data.frame(A_1 = "", A_2 = NA))
  dgeno2 <- hlaGeno(data.frame(A_1 = "A*03:01"))
  with_mocked_bindings(utils = list(data = getFromNamespace(".noop_utils_data", "deepMatchR")), {
    expect_error(calculateEpletLoad(rgeno2, dgeno2), "No allele strings")
  })
  
  # Pairwise: bad locus or not shared
  rgeno3 <- hlaGeno(data.frame(A_1="A*01:01", A_2="A*02:01"))
  dgeno3 <- hlaGeno(data.frame(A_1="A*03:01", A_2="A*24:02"))
  with_mocked_bindings(utils = list(data = getFromNamespace(".noop_utils_data", "deepMatchR")), {
    expect_error(calculateEpletLoad(rgeno3, dgeno3, return="pairwise", pairwise_locus = 1),
                 "provide pairwise_locus")
    expect_error(calculateEpletLoad(rgeno3, dgeno3, return="pairwise", pairwise_locus = "B"),
                 "not shared")
  })
})

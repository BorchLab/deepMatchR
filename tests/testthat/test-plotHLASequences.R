# tests/testthat/test-plotHLASequences.R

# Test data: mock sequences for testing
mock_seq_db <- list(
  "A*01:01" = "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEP",
  "A*02:01" = "MAVMAPRTLVLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQTDRVDLGTLRGYYNQSEAGSHIIQRMYGCDVGSDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQITKRKWEAARVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDPPKTHMTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEP",
  "A*03:01" = "MAVMAPRTLLLLLSGALALTQTWAGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQKMEPRAPWIEQEGPEYWDQETRNMKAHSQTDRANLGTLRGYYNQSEAGSHTLQSMYGCDVGPDGRFLRGYRQFAYDGKDYIALNEDLSSWTAADTAAQITQRKWEAARVAEQLRAYLEGECVEWLRRYLENGKDKLERADPPKTHVTHHPVSDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPKPLTLRWEP",
  "B*07:02" = "MRVTAPRTLLLLLWGAVALTETWAGSHSMRYFYTSVSRPGRGEPRFISVGYVDDTQFVRFDSDAASPRGEPRAPWVEQEGPEYWDRETQKYKRQAQADRVNLRKLRGYYNQSEDGSHTLQWMCGCDLGPDGRLLRGYDQYAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAARVAEQLRAYLEGLCVEWLRRYLENGKDTLERADPPKTHVTHHPVSDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWASVVVPSGQEQRYTCHMQHEGLPKPLTLRWEP",
  "B*08:01" = "MRVMAPRTLLLLLWGAVALTETWAGSHSMRYFYTAVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASPRGEPRAPWVEQEGPEYWDRNTQIFKTNTQTYRESLRNLRGYYNQSEAGSHTLQWMYGCDLGPDGRLLRGHNQYAYDGKDYIALNEDLRSWTAADKAAQITQRKWEAARVAEQLRAYLEGECVEWLRRYLENGKETLQRADPPKTHVTHHPISDHEATLRCWALGFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWASVVVPSGQEQRYTCHVQHEGLPKPLTLRWEP",
  "DRB1*03:01" = "MVCLKLPGGSCMTALTVTLMVLSSPLALAGDTRPRFLWQLKFECHFFNGTERVRLLERCIYNQEESVRFDSDVGEYRAVTELGRPDAEYWNSQKDLLEQRRAAVDTYCRHNYGVGESFTVQRRVHPEVTVYPAKTQPLQHHNLLVCSVSGFYPGSIEVRWFRNGQEEKTGVVSTGLIHNGDWTFQTLVMLETVPRSGEVYTCQVEHPSVTSPLTVEWRARSESAQSK",
  "DRB1*04:01" = "MVCLKLPGGSCMTALTVTLMVLSSPLALAGDTRPRFLEYSTSECHFFNGTERVRFLDRYFYHQEEYVRFDSDVGEFRAVTELGRPDAEYWNSQKDILEQARAAVDTYCRHNYGESFTVQRRVHPKVTVYPSKTQPLQHHNLLVCSVNGFYPGSIEVRWFRNGQEEKAGVVSTGLIQNGDWTFQTLVMLETVPRSGEVYTCQVEHPSLTSPLTVEWRARSESAQSK",
  "DRB1*07:01" = "MVCLKLPGGSCMTALTVTLMVLSSPLALAGDTRPRFLEQVKHECHFFNGTERVRLLERCIYNQEESVRFDSDVGEYRAVTELGRPDAEYWNSQKDLLEQRRAAVDTYCRHNYGVGESFTVQRRVHPEVTVYPAKTQPLQHHNLLVCSVSGFYPGSIEVRWFRNGQEEKTGVVSTGLIHNGDWTFQTLVMLETVPRSGEVYTCQVEHPSVTSPLTVEWRARSESAQSK"
)

# Create mock function for getAlleleSequence
mock_getAlleleSequence <- function(allele_name, type = "PROT", use_cache = TRUE, cache_dir = NULL) {
  if (type == "NUC") {
    # Return simple nucleotide sequences for testing
    nuc_seqs <- list(
      "A*01:01" = "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTACACCTCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAG",
      "A*02:01" = "ATGGCCGTCATGGCGCCCCGAACCCTCGTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTACACCTCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAG"
    )
    if (allele_name %in% names(nuc_seqs)) {
      return(nuc_seqs[[allele_name]])
    }
    stop(sprintf("Allele '%s' not found in mock database for type 'NUC'.", allele_name))
  }

  if (allele_name %in% names(mock_seq_db)) {
    return(mock_seq_db[[allele_name]])
  }
  stop(sprintf("Allele '%s' not found in mock IMGT/HLA database for type 'PROT'.", allele_name))
}

test_that("plotHLASequences works with two protein sequences", {
  testthat::with_mocked_bindings(
    getAlleleSequence = mock_getAlleleSequence,
    {
      result <- plotHLASequences(
        alleles = c("A*01:01", "A*02:01"),
        seq_type = "protein",
        verbose = FALSE
      )

      expect_s3_class(result, "hla_sequence_plot")
      expect_true("plot" %in% names(result))
      expect_true("mismatch_summary" %in% names(result))
      expect_true("mismatch_details" %in% names(result))
      expect_true("sequences" %in% names(result))

      # Check mismatch summary has correct structure
      expect_s3_class(result$mismatch_summary, "data.frame")
      expect_equal(nrow(result$mismatch_summary), 1)  # One comparison
      expect_true(all(c("comparison", "total_mismatches", "filtered_mismatches",
                        "comparable_positions", "percent_identity") %in%
                        names(result$mismatch_summary)))
    }
  )
})

test_that("plotHLASequences works with three or more sequences", {
  testthat::with_mocked_bindings(
    getAlleleSequence = mock_getAlleleSequence,
    {
      result <- plotHLASequences(
        alleles = c("A*01:01", "A*02:01", "A*03:01"),
        seq_type = "protein",
        verbose = FALSE
      )

      expect_s3_class(result, "hla_sequence_plot")
      expect_equal(length(result$sequences), 3)

      # Should have 2 comparisons (all vs reference)
      expect_equal(nrow(result$mismatch_summary), 2)
    }
  )
})

test_that("plotHLASequences works with custom sequences", {
  custom_seqs <- list(
    "Seq1" = "YFAMYGEKVAHTHVDTLYVRYHY",
    "Seq2" = "YFDMYGEKVAHTHVDTLYVRFHY",
    "Seq3" = "YFAMYGEKVAHTHVDTLYVRFHY"
  )

  result <- plotHLASequences(
    alleles = names(custom_seqs),
    sequences = custom_seqs,
    seq_type = "protein",
    show_exons = FALSE,
    verbose = FALSE
  )

  expect_s3_class(result, "hla_sequence_plot")
  expect_equal(length(result$sequences), 3)
  expect_equal(names(result$sequences), names(custom_seqs))
})

test_that("plotHLASequences respects filter_charge parameter", {
  custom_seqs <- list(
    "Ref" = "AAAS",  # A=nonpolar, A=nonpolar, A=nonpolar, S=polar
    "Alt" = "ADST"   # A=nonpolar, D=negative, S=polar, T=polar
  )
  # Mismatches: A->D (charge change), A->S (no charge change), S->T (no charge change)

  result_all <- plotHLASequences(
    alleles = names(custom_seqs),
    sequences = custom_seqs,
    seq_type = "protein",
    show_exons = FALSE,
    filter_charge = NULL,
    verbose = FALSE
  )

  result_charge <- plotHLASequences(
    alleles = names(custom_seqs),
    sequences = custom_seqs,
    seq_type = "protein",
    show_exons = FALSE,
    filter_charge = TRUE,
    verbose = FALSE
  )

  # Without filter: 3 mismatches
  expect_equal(result_all$mismatch_summary$total_mismatches[1], 3)

  # With charge filter: only 1 mismatch (A->D changes charge)
  expect_equal(result_charge$mismatch_summary$filtered_mismatches[1], 1)
})

test_that("plotHLASequences respects filter_polarity parameter", {
  custom_seqs <- list(
    "Ref" = "AAAS",
    "Alt" = "ADST"
  )
  # A->D: polarity change (nonpolar->polar)
  # A->S: polarity change (nonpolar->polar)
  # S->T: no polarity change (both polar)

  result_polarity <- plotHLASequences(
    alleles = names(custom_seqs),
    sequences = custom_seqs,
    seq_type = "protein",
    show_exons = FALSE,
    filter_polarity = TRUE,
    verbose = FALSE
  )

  # With polarity filter: 2 mismatches change polarity
  expect_equal(result_polarity$mismatch_summary$filtered_mismatches[1], 2)
})

test_that("plotHLASequences handles reference_idx parameter", {
  custom_seqs <- list(
    "Seq1" = "AAAA",
    "Seq2" = "AAAD",
    "Seq3" = "AADA"
  )

  # Reference is Seq1 (default)
  result1 <- plotHLASequences(
    alleles = names(custom_seqs),
    sequences = custom_seqs,
    seq_type = "protein",
    show_exons = FALSE,
    reference_idx = 1,
    verbose = FALSE
  )

  # Reference is Seq2
  result2 <- plotHLASequences(
    alleles = names(custom_seqs),
    sequences = custom_seqs,
    seq_type = "protein",
    show_exons = FALSE,
    reference_idx = 2,
    verbose = FALSE
  )

  # Comparisons should be different
  expect_true(grepl("Seq1 vs", result1$mismatch_summary$comparison[1]))
  expect_true(grepl("Seq2 vs", result2$mismatch_summary$comparison[1]))
})

test_that("plotHLASequences includes exon boundaries for known loci", {
  testthat::with_mocked_bindings(
    getAlleleSequence = mock_getAlleleSequence,
    {
      result <- plotHLASequences(
        alleles = c("DRB1*03:01", "DRB1*04:01"),
        seq_type = "protein",
        show_exons = TRUE,
        verbose = FALSE
      )

      # DRB1 should have known exon boundaries
      expect_false(is.null(result$exon_boundaries))
      expect_false(is.null(result$exon_descriptions))
      expect_true(length(result$exon_boundaries) > 0)
    }
  )
})

test_that("plotHLASequences returns correct plot components", {
  custom_seqs <- list(
    "Seq1" = "YFAMYGEKVAHTHVDTLYVRYHY",
    "Seq2" = "YFDMYGEKVAHTHVDTLYVRFHY"
  )

  result <- plotHLASequences(
    alleles = names(custom_seqs),
    sequences = custom_seqs,
    seq_type = "protein",
    show_exons = FALSE,
    verbose = FALSE
  )

  # Should have main plot
  expect_true("main" %in% names(result$plots))

  # Should have mismatch plot (since there are mismatches)
  expect_true("mismatch" %in% names(result$plots))

  # Should have legend for protein with properties color scheme
  expect_true("legend" %in% names(result$plots))
})

test_that("plotHLASequences handles focus_on_mismatches parameter", {
  custom_seqs <- list(
    "Seq1" = "YFAMYGEKVAHTHVDTLYVRYHY",
    "Seq2" = "YFDMYGEKVAHTHVDTLYVRFHY"
  )

  result <- plotHLASequences(
    alleles = names(custom_seqs),
    sequences = custom_seqs,
    seq_type = "protein",
    show_exons = FALSE,
    focus_on_mismatches = TRUE,
    context_window = 5,
    verbose = FALSE
  )

  # Should have zoom plot when focusing on mismatches
  expect_true("zoom" %in% names(result$plots))
})

test_that("plotHLASequences errors with less than 2 sequences", {
  expect_error(
    plotHLASequences(
      alleles = c("A*01:01"),
      seq_type = "protein",
      verbose = FALSE
    ),
    "At least 2 alleles"
  )
})

test_that("plotHLASequences errors with invalid reference_idx", {
  custom_seqs <- list(
    "Seq1" = "AAAA",
    "Seq2" = "AAAD"
  )

  expect_error(
    plotHLASequences(
      alleles = names(custom_seqs),
      sequences = custom_seqs,
      seq_type = "protein",
      reference_idx = 5,
      verbose = FALSE
    ),
    "reference_idx must be between"
  )
})

test_that("plotHLASequences handles identical sequences", {
  custom_seqs <- list(
    "Seq1" = "AAAA",
    "Seq2" = "AAAA"
  )

  result <- plotHLASequences(
    alleles = names(custom_seqs),
    sequences = custom_seqs,
    seq_type = "protein",
    show_exons = FALSE,
    verbose = FALSE
  )

  expect_equal(result$mismatch_summary$total_mismatches[1], 0)
  expect_equal(result$mismatch_summary$percent_identity[1], 100)
})

test_that("plotHLASequences returns alignment_df with correct structure", {
  custom_seqs <- list(
    "Seq1" = "AAAA",
    "Seq2" = "AAAD"
  )

  result <- plotHLASequences(
    alleles = names(custom_seqs),
    sequences = custom_seqs,
    seq_type = "protein",
    show_exons = FALSE,
    verbose = FALSE
  )

  expect_s3_class(result$alignment_df, "data.frame")
  expect_true(all(c("Position", "Sequence", "Residue") %in% names(result$alignment_df)))

  # Should have rows for each sequence at each position
  expect_equal(nrow(result$alignment_df), 4 * 2)  # 4 positions x 2 sequences
})

test_that("plotHLASequences handles different color schemes", {
  custom_seqs <- list(
    "Seq1" = "YFAMY",
    "Seq2" = "YFDMY"
  )

  for (scheme in c("properties", "hydropathy", "classic")) {
    result <- plotHLASequences(
      alleles = names(custom_seqs),
      sequences = custom_seqs,
      seq_type = "protein",
      show_exons = FALSE,
      color_scheme = scheme,
      verbose = FALSE
    )

    expect_s3_class(result, "hla_sequence_plot")
  }
})

test_that("print.hla_sequence_plot works", {
  custom_seqs <- list(
    "Seq1" = "AAAA",
    "Seq2" = "AAAD"
  )

  result <- plotHLASequences(
    alleles = names(custom_seqs),
    sequences = custom_seqs,
    seq_type = "protein",
    show_exons = FALSE,
    verbose = FALSE
  )

  # Print method should return invisibly
  expect_invisible(print(result))
})

test_that("plotHLASequences works with nucleotide sequences", {
  testthat::with_mocked_bindings(
    getAlleleSequence = mock_getAlleleSequence,
    {
      result <- plotHLASequences(
        alleles = c("A*01:01", "A*02:01"),
        seq_type = "nucleotide",
        verbose = FALSE
      )

      expect_s3_class(result, "hla_sequence_plot")
      expect_true("plot" %in% names(result))
    }
  )
})

test_that("internal helper .inferLocus works correctly", {
  # Access internal function
  inferLocus <- deepMatchR:::.inferLocus

 expect_equal(inferLocus("A*01:01"), "A")
  expect_equal(inferLocus("DRB1*03:01"), "DRB1")
  expect_equal(inferLocus("DQB1*02:01"), "DQB1")
  expect_equal(inferLocus("B*07:02"), "B")
})

test_that("mismatch_details contains domain information when available", {
  testthat::with_mocked_bindings(
    getAlleleSequence = mock_getAlleleSequence,
    {
      result <- plotHLASequences(
        alleles = c("DRB1*03:01", "DRB1*04:01"),
        seq_type = "protein",
        show_exons = TRUE,
        verbose = FALSE
      )

      # Check that domain column exists in mismatch details
      detail_df <- result$mismatch_details[[1]]
      expect_true("domain" %in% names(detail_df))
    }
  )
})

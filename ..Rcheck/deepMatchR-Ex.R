pkgname <- "deepMatchR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('deepMatchR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("compareHlaSequences")
### * compareHlaSequences

flush(stderr()); flush(stdout())

### Name: compareHlaSequences
### Title: Compare two HLA protein sequences and identify polymorphisms
### Aliases: compareHlaSequences

### ** Examples

seq1 <- "YFAMYGEKVAHTHVDTLYVRYHY"
seq2 <- "YFDMYGEKVAHTHVDTLYVRYHY"
compareHlaSequences(seq1, seq2)



cleanEx()
nameEx("plotEplets")
### * plotEplets

flush(stderr()); flush(stdout())

### Name: plotEplets
### Title: Plot Eplet Results from SPI Assay
### Aliases: plotEplets

### ** Examples

# Using a data frame:
plotEplets(deepMatchR_example[[1]],
           cutoff = 2000,
           evidence_level = c("A1", "A2", "B"),
           percPos_filter = 0.4,
           plot_type = "treemap")



cleanEx()
nameEx("quantifyMismatch")
### * quantifyMismatch

flush(stderr()); flush(stdout())

### Name: quantifyMismatch
### Title: Quantify Amino Acid Mismatches Between Two Sequences
### Aliases: quantifyMismatch

### ** Examples

seq1 <- "YFAMYGEKVAHTHVDTLYVRYHY"
seq2 <- "YFDMYGEKVAHTHVDTLYVRFHY"
quantifyMismatch(seq1, seq2)
#> [1] 2




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

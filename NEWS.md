# deepMatchR 0.99.0

## Initial Release

* Core functionality for HLA matching and immunogenicity assessment
* Sequence retrieval with `getAlleleSequence()` and `batchGetSequences()`
* Mismatch quantification with `quantifyMismatch()` including biophysical filters
* Mismatch load calculation with `calculateMismatchLoad()`
* Eplet analysis with `quantifyEpletMismatch()` and `calculateEpletLoad()`
* Peptide binding prediction support (MHCnuggets, NetMHCpan, PWM backends)
* HLA genotype handling with `hlaGeno()`
* Added `plotHLASequences()` function for visualizing and comparing two or more
  HLA allele sequences (nucleotide or protein)
* Visualization functions for antibodies and eplets

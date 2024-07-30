#Evaluate the MHC-I or MHC-II mismatch score

#TODO Evaluate sequences for mismatch
#TODO Extract the flanking sequences from reference
#TODO Feed sequences into mhcnuggets
#TODO Develop harmonic system



#Make all the combinations of peptides with the polymorphic regions
.extract.combinations <- function(input_string, start_pos, length) {
  if(start_pos > nchar(input_string) || start_pos < 1) {
    stop("Invalid start position")
  }
  if(length < 1 || length > nchar(input_string)) {
    stop("Invalid length")
  }
  
  substrings <- c()
  start_limit <- max(1, start_pos - length + 1)
  
  for(i in start_limit:(nchar(input_string) - length + 1)) {
    substrings <- c(substrings, substr(input_string, i, i + length - 1))
  }
  
  return(substrings)
}

# Example usage:
input_string <- "AANAAA"
start_pos <- 3
length <- 4
extract_combinations(input_string, start_pos, length)

#add correlation feature to pin files during entrapment runs

#packages
library(tidyverse)

#get arguments
args = commandArgs(trailingOnly = TRUE)

#function
get_correlation_feature = function(file, p, file_path) {
  df = read_delim(file)
  df$corr_feature = 0
  df$corr_feature[df$Label == 1] = sample(c(0, 1), size = sum(df$Label == 1), replace = TRUE, prob = c(1 - p, p))
  write_delim(df, file_path, delim = '\t', quote = "none")
}

#run
get_correlation_feature(as.character(args[1]), as.numeric(args[2]), as.character(args[3]))


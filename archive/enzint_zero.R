#set the args
args = commandArgs(trailingOnly=TRUE)
enzint_zero = function(file_name, name) {
  print('success')
  
  #libraries
  library(tidyverse)
  options(scipen = 999)
  
  print('deleting enzInt')
  
  
  file = read_delim(file_name, na = character())
  file$enzInt = 0
  write.table(file, paste(name, sep = ''), row.names = F, sep = '\t', quote = F)
  
}

enzint_zero(args[1], args[2])
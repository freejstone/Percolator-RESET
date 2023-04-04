#set the args
args = commandArgs(trailingOnly=TRUE)
enzint_zero = function(file_name, file_name_2, dir, count) {
  print('success')
  
  #libraries
  library(tidyverse)
  options(scipen = 999)
  
  print('doing splitting')
  
  
  file = read_delim(file_name, na = character())
  file$enzInt = 0
  file_2 = read_delim(file_name_2, na = character())
  file_2 = file_2[file_2$`decoy index` == count, ]
  unique_id = paste(file_2$scan, file_2$sequence)
  Peptide = file$Peptide
  Peptide = gsub(".*[.]([^.]+)[.].*", "\\1", Peptide)
  unique_id_file = paste(file$ScanNr, Peptide)
  file = file[unique_id_file %in% unique_id, ]
  write.table(file, paste(dir, '-', count, '.tide-search.decoy.pin', sep = ''), row.names = F, sep = '\t', quote = F)
  
}

enzint_zero(args[1], args[2], args[3], args[4])
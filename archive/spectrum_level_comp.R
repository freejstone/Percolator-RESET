#set the args
args = commandArgs(trailingOnly=TRUE)
spectrum_level_comp = function(target_pin, decoy_pin, name) {
  print('success')
  
  #libraries
  library(tidyverse)
  options(scipen = 999)
  
  print('doing spectrum-level competition')
  
  target = read_delim(target_pin, na = character())
  decoy = read_delim(decoy_pin, na = character())
  target_decoy = rbind(target, decoy)
  target_decoy = target_decoy[sample(nrow(target_decoy)), ]
  target_decoy = target_decoy[order(-target_decoy$TailorScore), ]
  target_decoy = target_decoy[!duplicated(target_decoy$ScanNr), ]
  
  write.table(target_decoy, paste(name, '.tide-search.pin', sep = ''), row.names = F, sep = '\t', quote = F)
  
}

spectrum_level_comp(args[1], args[2], args[3])
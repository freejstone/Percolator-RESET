#the following produces warning, but they can be safely ignored

options(scipen = 999) #to read/write floats to the exact precision
library(tidyverse)
for (d in 0:99) {
  search_output = read_delim(paste('tide-search/full/open_top5_full_', d, '.tide-search.pin', sep = ''), delim = '\t', show_col_types = FALSE)
  search_output$enzInt = 0
  rank = unlist(lapply(search_output$SpecId, function(x) substr(x, nchar(x), nchar(x))))
  search_output = search_output[rank == "1", ]
  peptide_list = read.delim(paste('tide-index/index-full-', d, '/tide-index.peptides.txt', sep = ''))
  low_complex = peptide_list$target[peptide_list$decoy.s. == '']
  peptide_no_flank = unlist(lapply(search_output$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
  print(sum(peptide_no_flank %in% low_complex))
  search_output = search_output[!(peptide_no_flank %in% low_complex), ]
  write.table(search_output, paste('tide-search/full/open_top1_enzint_zero_full_', d, '.tide-search.pin', sep = ''), row.names = F, sep = '\t', quote = F)
}

####################################################################################

for (d in 0:99) {
  search_output = read_delim(paste('tide-search/full/dcy_full_', d, '.tide-search.pin', sep = ''), delim = '\t', show_col_types = FALSE)
  search_output$enzInt = 0
  rank = unlist(lapply(search_output$SpecId, function(x) substr(x, nchar(x), nchar(x))))
  search_output = search_output[rank == "1", ]
  peptide_list = read.delim(paste('tide-index/index-full-', d, '/tide-index.peptides.txt', sep = ''))
  low_complex = peptide_list$target[peptide_list$decoy.s. == '']
  peptide_no_flank = unlist(lapply(search_output$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
  search_output = search_output[!(peptide_no_flank %in% low_complex), ]
  write.table(search_output, paste('tide-search/full/dcy_enzint_zero_full_', d, '.tide-search.pin', sep = ''), row.names = F, sep = '\t', quote = F)
}

####################################################################################

for (d in 0:99) {
  search_output = read_delim(paste('tide-search/half/open_top5_half_', d, '.tide-search.pin', sep = ''), delim = '\t', show_col_types = FALSE)
  search_output$enzInt = 0
  rank = unlist(lapply(search_output$SpecId, function(x) substr(x, nchar(x), nchar(x))))
  search_output = search_output[rank == "1", ]
  peptide_list = read.delim(paste('tide-index/index-half-', d, '/tide-index.peptides.txt', sep = ''))
  low_complex = peptide_list$target[peptide_list$decoy.s. == '']
  peptide_no_flank = unlist(lapply(search_output$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
  search_output = search_output[!(peptide_no_flank %in% low_complex), ]
  write.table(search_output, paste('tide-search/half/open_top1_enzint_zero_half_', d, '.tide-search.pin', sep = ''), row.names = F, sep = '\t', quote = F)
}

####################################################################################

for (d in 0:99) {
  search_output = read_delim(paste('tide-search/half/dcy_half_', d, '.tide-search.pin', sep = ''), delim = '\t', show_col_types = FALSE)
  search_output$enzInt = 0
  rank = unlist(lapply(search_output$SpecId, function(x) substr(x, nchar(x), nchar(x))))
  search_output = search_output[rank == "1", ]
  peptide_list = read.delim(paste('tide-index/index-half-', d, '/tide-index.peptides.txt', sep = ''))
  low_complex = peptide_list$target[peptide_list$decoy.s. == '']
  peptide_no_flank = unlist(lapply(search_output$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
  search_output = search_output[!(peptide_no_flank %in% low_complex), ]
  write.table(search_output, paste('tide-search/half/dcy_enzint_zero_half_', d, '.tide-search.pin', sep = ''), row.names = F, sep = '\t', quote = F)
}

####################################################################################

for (d in 0:99) {
  search_output = read_delim(paste('tide-search/025/open_top5_025_', d, '.tide-search.pin', sep = ''), delim = '\t', show_col_types = FALSE)
  search_output$enzInt = 0
  rank = unlist(lapply(search_output$SpecId, function(x) substr(x, nchar(x), nchar(x))))
  search_output = search_output[rank == "1", ]
  peptide_list = read.delim(paste('tide-index/index-025-', d, '/tide-index.peptides.txt', sep = ''))
  low_complex = peptide_list$target[peptide_list$decoy.s. == '']
  peptide_no_flank = unlist(lapply(search_output$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
  search_output = search_output[!(peptide_no_flank %in% low_complex), ]
  write.table(search_output, paste('tide-search/025/open_top1_enzint_zero_025_', d, '.tide-search.pin', sep = ''), row.names = F, sep = '\t', quote = F)
}

####################################################################################

for (d in 0:99) {
  search_output = read_delim(paste('tide-search/025/dcy_025_', d, '.tide-search.pin', sep = ''), delim = '\t', show_col_types = FALSE)
  search_output$enzInt = 0
  rank = unlist(lapply(search_output$SpecId, function(x) substr(x, nchar(x), nchar(x))))
  search_output = search_output[rank == "1", ]
  peptide_list = read.delim(paste('tide-index/index-025-', d, '/tide-index.peptides.txt', sep = ''))
  low_complex = peptide_list$target[peptide_list$decoy.s. == '']
  peptide_no_flank = unlist(lapply(search_output$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
  search_output = search_output[!(peptide_no_flank %in% low_complex), ]
  write.table(search_output, paste('tide-search/025/dcy_enzint_zero_025_', d, '.tide-search.pin', sep = ''), row.names = F, sep = '\t', quote = F)
}

####################################################################################

for (d in 0:99) {
  search_output = read_delim(paste('tide-search/075/open_top5_075_', d, '.tide-search.pin', sep = ''), delim = '\t', show_col_types = FALSE)
  search_output$enzInt = 0
  rank = unlist(lapply(search_output$SpecId, function(x) substr(x, nchar(x), nchar(x))))
  search_output = search_output[rank == "1", ]
  peptide_list = read.delim(paste('tide-index/index-075-', d, '/tide-index.peptides.txt', sep = ''))
  low_complex = peptide_list$target[peptide_list$decoy.s. == '']
  peptide_no_flank = unlist(lapply(search_output$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
  search_output = search_output[!(peptide_no_flank %in% low_complex), ]
  write.table(search_output, paste('tide-search/075/open_top1_enzint_zero_075_', d, '.tide-search.pin', sep = ''), row.names = F, sep = '\t', quote = F)
}

####################################################################################

for (d in 0:99) {
  search_output = read_delim(paste('tide-search/075/dcy_075_', d, '.tide-search.pin', sep = ''), delim = '\t', show_col_types = FALSE)
  search_output$enzInt = 0
  rank = unlist(lapply(search_output$SpecId, function(x) substr(x, nchar(x), nchar(x))))
  search_output = search_output[rank == "1", ]
  peptide_list = read.delim(paste('tide-index/index-075-', d, '/tide-index.peptides.txt', sep = ''))
  low_complex = peptide_list$target[peptide_list$decoy.s. == '']
  peptide_no_flank = unlist(lapply(search_output$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
  search_output = search_output[!(peptide_no_flank %in% low_complex), ]
  write.table(search_output, paste('tide-search/075/dcy_enzint_zero_075_', d, '.tide-search.pin', sep = ''), row.names = F, sep = '\t', quote = F)
}

####################################################################################
# 
# 
# 
# 
# 
# 
# 
# 
# options(scipen = 999)
# library(tidyverse)
# for (d in 0:99) {
#   search_output = read_delim(paste('full/open_top5_sep_full_', d, '.tide-search.target.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('full/open_top1_sep_enzint_zero_full_', d, '.tide-search.target.pin', sep = ''), row.names = F, sep = '\t', quote = F)
#   
#   search_output = read_delim(paste('full/open_top5_sep_full_', d, '.tide-search.decoy.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('full/open_top1_sep_enzint_zero_full_', d, '.tide-search.decoy.pin', sep = ''), row.names = F, sep = '\t', quote = F)
# }
# 
# ####################################################################################
# 
# options(scipen = 999)
# library(tidyverse)
# for (d in 0:99) {
#   search_output = read_delim(paste('full/dcy_sep_full_', d, '.tide-search.target.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('full/dcy_sep_enzint_zero_full_', d, '.tide-search.target.pin', sep = ''), row.names = F, sep = '\t', quote = F)
#   
#   search_output = read_delim(paste('full/dcy_sep_full_', d, '.tide-search.decoy.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('full/dcy_sep_enzint_zero_full_', d, '.tide-search.decoy.pin', sep = ''), row.names = F, sep = '\t', quote = F)
# }
# ####################################################################################
# 
# 
# 
# options(scipen = 999)
# library(tidyverse)
# for (d in 0:99) {
#   search_output = read_delim(paste('075/open_top5_sep_075_', d, '.tide-search.target.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('075/open_top1_sep_enzint_zero_075_', d, '.tide-search.target.pin', sep = ''), row.names = F, sep = '\t', quote = F)
#   
#   search_output = read_delim(paste('075/open_top5_sep_075_', d, '.tide-search.decoy.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('075/open_top1_sep_enzint_zero_075_', d, '.tide-search.decoy.pin', sep = ''), row.names = F, sep = '\t', quote = F)
# }
# 
# ####################################################################################
# 
# options(scipen = 999)
# library(tidyverse)
# for (d in 0:99) {
#   search_output = read_delim(paste('075/dcy_sep_075_', d, '.tide-search.target.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('075/dcy_sep_enzint_zero_075_', d, '.tide-search.target.pin', sep = ''), row.names = F, sep = '\t', quote = F)
#   
#   search_output = read_delim(paste('075/dcy_sep_075_', d, '.tide-search.decoy.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('075/dcy_sep_enzint_zero_075_', d, '.tide-search.decoy.pin', sep = ''), row.names = F, sep = '\t', quote = F)
# }


#########################


# 
# 
# 
# options(scipen = 999)
# library(tidyverse)
# for (d in 0:99) {
#   search_output = read_delim(paste('half/open_top5_sep_half_', d, '.tide-search.target.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('half/open_top1_sep_enzint_zero_half_', d, '.tide-search.target.pin', sep = ''), row.names = F, sep = '\t', quote = F)
#   
#   search_output = read_delim(paste('half/open_top5_sep_half_', d, '.tide-search.decoy.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('half/open_top1_sep_enzint_zero_half_', d, '.tide-search.decoy.pin', sep = ''), row.names = F, sep = '\t', quote = F)
# }
# 
# ####################################################################################
# 
# options(scipen = 999)
# library(tidyverse)
# for (d in 0:99) {
#   search_output = read_delim(paste('half/dcy_sep_half_', d, '.tide-search.target.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('half/dcy_sep_enzint_zero_half_', d, '.tide-search.target.pin', sep = ''), row.names = F, sep = '\t', quote = F)
#   
#   search_output = read_delim(paste('half/dcy_sep_half_', d, '.tide-search.decoy.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('half/dcy_sep_enzint_zero_half_', d, '.tide-search.decoy.pin', sep = ''), row.names = F, sep = '\t', quote = F)
# }
# 
# 
# ##########
# 
# 
# 
# 
# 
# options(scipen = 999)
# library(tidyverse)
# for (d in 0:99) {
#   search_output = read_delim(paste('025/open_top5_sep_025_', d, '.tide-search.target.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('025/open_top1_sep_enzint_zero_025_', d, '.tide-search.target.pin', sep = ''), row.names = F, sep = '\t', quote = F)
#   
#   search_output = read_delim(paste('025/open_top5_sep_025_', d, '.tide-search.decoy.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('025/open_top1_sep_enzint_zero_025_', d, '.tide-search.decoy.pin', sep = ''), row.names = F, sep = '\t', quote = F)
# }
# 
# ####################################################################################
# 
# options(scipen = 999)
# library(tidyverse)
# for (d in 0:99) {
#   search_output = read_delim(paste('025/dcy_sep_025_', d, '.tide-search.target.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('025/dcy_sep_enzint_zero_025_', d, '.tide-search.target.pin', sep = ''), row.names = F, sep = '\t', quote = F)
#   
#   search_output = read_delim(paste('025/dcy_sep_025_', d, '.tide-search.decoy.pin', sep = ''), delim = '\t', show_col_types = FALSE)
#   search_output$enzInt = 0
#   unique_ids = gsub('.{2}$', '', search_output$SpecId)
#   unique_ids = gsub('target', '', unique_ids)
#   unique_ids = gsub('decoy', '', unique_ids)
#   search_output$unique_ids = unique_ids
#   search_output = search_output[order(search_output$unique_ids, -search_output$XCorr), ] #preserves the order of ties within the search file
#   search_output = search_output %>% group_by(unique_ids) %>% mutate(rank = order(unique_ids)) %>% ungroup()
#   search_output = search_output %>% filter(rank %in% 1:1)
#   search_output = search_output %>% select(-unique_ids, -rank)
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'decoy_ENTR'
#   search_output$Proteins[grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'decoy_ISB18'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ENTR', search_output$Proteins)] = 'target_ENTR'
#   search_output$Proteins[!grepl('decoy', search_output$Proteins) & grepl('ISB18', search_output$Proteins)] = 'target_ISB18'
#   
#   write.table(search_output, paste('025/dcy_sep_enzint_zero_025_', d, '.tide-search.decoy.pin', sep = ''), row.names = F, sep = '\t', quote = F)
# }
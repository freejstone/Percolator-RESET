#combine and create fasta file from peptide lists
library(tidyverse)
options(scipen = 999)
peptide_list_ISB18 = read.delim('./ISB18/tide-index.peptides.txt')
peptide_list_yeast = read.delim('./castor/tide-index.peptides.txt')

peptide_list_ISB18$pr = paste(peptide_list_ISB18$pr,'_ISB18', sep = '')
peptide_list_yeast$pr = paste(peptide_list_yeast$pr,'_ENTR', sep = '')

targets = c(peptide_list_ISB18$target, peptide_list_yeast$target)
proteins = c(peptide_list_ISB18$pr, peptide_list_yeast$pr)

peptide_list = data.frame(targets, proteins)
names(peptide_list) = c('targets', 'proteins')

peptide_list = peptide_list %>% distinct %>% 
  arrange(targets, proteins) %>% 
  group_by(targets) %>%
  summarise(proteins = str_c(proteins, collapse=","))

peptide_list = peptide_list[!(grepl('ISB18', peptide_list$proteins) & grepl('ENTR', peptide_list$proteins)), ]

peptide_list$target_L = gsub('I', 'L', peptide_list$targets)
ISB18_target_L = peptide_list$target_L[grepl('ISB18', peptide_list$proteins)]
peptide_list_sub = peptide_list %>% 
  filter(target_L %in% ISB18_target_L) %>% group_by(target_L) %>% mutate(check = if (n() == 1) { FALSE } else { any(grepl('ISB18', proteins)) & any(grepl('ENTR', proteins)) } ) %>% 
  ungroup() %>% filter(check == TRUE) %>% filter(grepl('ENTR', proteins))
peptide_list = peptide_list %>% filter(!(targets %in% peptide_list_sub$targets))

peptide_list$proteins = paste(">", peptide_list$proteins, sep = "") #898 ISB18 peptides, 571298 ENTR peptides

peptide_list = peptide_list[, 2:1]

vec_pep_list = as.vector(t(as.matrix(peptide_list)))
file = file("combine_full.fasta")
writeLines(vec_pep_list, file)


########################################################################################################################

#combine and create fasta file from peptide lists
library(tidyverse)
options(scipen = 999)
peptide_list_ISB18 = read.delim('./ISB18/tide-index.peptides.txt')
peptide_list_yeast = read.delim('./castor/tide-index.peptides.txt')

peptide_list_ISB18$pr = paste(peptide_list_ISB18$pr,'_ISB18', sep = '')

peptide_list_yeast$pr = paste(peptide_list_yeast$pr,'_ENTR', sep = '')

targets = c(peptide_list_ISB18$target, peptide_list_yeast$target)
proteins = c(peptide_list_ISB18$pr, peptide_list_yeast$pr)

peptide_list = data.frame(targets, proteins)
names(peptide_list) = c('targets', 'proteins')

peptide_list = peptide_list %>% distinct %>% 
  arrange(targets, proteins) %>% 
  group_by(targets) %>%
  summarise(proteins = str_c(proteins, collapse=","))

#924 ISB18, 571333 ENTR

peptide_list = peptide_list[!(grepl('ISB18', peptide_list$proteins) & grepl('ENTR', peptide_list$proteins)), ]

peptide_list$target_L = gsub('I', 'L', peptide_list$targets)
ISB18_target_L = peptide_list$target_L[grepl('ISB18', peptide_list$proteins)]
peptide_list_sub = peptide_list %>% 
  filter(target_L %in% ISB18_target_L) %>% group_by(target_L) %>% mutate(check = if (n() == 1) { FALSE } else { any(grepl('ISB18', proteins)) & any(grepl('ENTR', proteins)) } ) %>% 
  ungroup() %>% filter(check == TRUE) %>% filter(grepl('ENTR', proteins))
peptide_list = peptide_list %>% filter(!(targets %in% peptide_list_sub$targets))

peptide_list$proteins = paste(">", peptide_list$proteins, sep = "") #898 ISB18 peptides, 571298 ENTR peptides

peptide_list = peptide_list[, 2:1]

ISB18_inds = which(grepl('ISB18', peptide_list$proteins))
set.seed(1)
ISB18_inds_remove = sample(ISB18_inds, size = length(ISB18_inds)/2)
peptide_list = peptide_list[-ISB18_inds_remove, ] #449 ISB18 peptides, 571298 ENTR peptides

vec_pep_list = as.vector(t(as.matrix(peptide_list)))
file = file("combine_half.fasta")
writeLines(vec_pep_list, file)
close(file)

########################################################################################################################

#combine and create fasta file from peptide lists
library(tidyverse)
options(scipen = 999)
peptide_list_ISB18 = read.delim('./ISB18/tide-index.peptides.txt')
peptide_list_yeast = read.delim('./castor/tide-index.peptides.txt')

peptide_list_ISB18$pr = paste(peptide_list_ISB18$pr,'_ISB18', sep = '')

peptide_list_yeast$pr = paste(peptide_list_yeast$pr,'_ENTR', sep = '')

targets = c(peptide_list_ISB18$target, peptide_list_yeast$target)
proteins = c(peptide_list_ISB18$pr, peptide_list_yeast$pr)

peptide_list = data.frame(targets, proteins)
names(peptide_list) = c('targets', 'proteins')

peptide_list = peptide_list %>% distinct %>% 
  arrange(targets, proteins) %>% 
  group_by(targets) %>%
  summarise(proteins = str_c(proteins, collapse=","))

#924 ISB18, 571333 ENTR

peptide_list = peptide_list[!(grepl('ISB18', peptide_list$proteins) & grepl('ENTR', peptide_list$proteins)), ]

peptide_list$target_L = gsub('I', 'L', peptide_list$targets)
ISB18_target_L = peptide_list$target_L[grepl('ISB18', peptide_list$proteins)]
peptide_list_sub = peptide_list %>% 
  filter(target_L %in% ISB18_target_L) %>% group_by(target_L) %>% mutate(check = if (n() == 1) { FALSE } else { any(grepl('ISB18', proteins)) & any(grepl('ENTR', proteins)) } ) %>% 
  ungroup() %>% filter(check == TRUE) %>% filter(grepl('ENTR', proteins))
peptide_list = peptide_list %>% filter(!(targets %in% peptide_list_sub$targets))

peptide_list$proteins = paste(">", peptide_list$proteins, sep = "") #898 ISB18 peptides, 571298 ENTR peptides

peptide_list = peptide_list[, 2:1]

ISB18_inds = which(grepl('ISB18', peptide_list$proteins))
set.seed(1)
ISB18_inds_remove = sample(ISB18_inds, size = floor(length(ISB18_inds)/4))
peptide_list = peptide_list[-ISB18_inds_remove, ] #674 ISB18 peptides, 571298 ENTR peptides

vec_pep_list = as.vector(t(as.matrix(peptide_list)))
file = file("combine_075.fasta")
writeLines(vec_pep_list, file)
close(file)
########################################################################################################################

#combine and create fasta file from peptide lists
library(tidyverse)
options(scipen = 999)
peptide_list_ISB18 = read.delim('./ISB18/tide-index.peptides.txt')
peptide_list_yeast = read.delim('./castor/tide-index.peptides.txt')

peptide_list_ISB18$pr = paste(peptide_list_ISB18$pr,'_ISB18', sep = '')

peptide_list_yeast$pr = paste(peptide_list_yeast$pr,'_ENTR', sep = '')

targets = c(peptide_list_ISB18$target, peptide_list_yeast$target)
proteins = c(peptide_list_ISB18$pr, peptide_list_yeast$pr)

peptide_list = data.frame(targets, proteins)
names(peptide_list) = c('targets', 'proteins')

peptide_list = peptide_list %>% distinct %>% 
  arrange(targets, proteins) %>% 
  group_by(targets) %>%
  summarise(proteins = str_c(proteins, collapse=","))

#924 ISB18, 571333 ENTR

peptide_list = peptide_list[!(grepl('ISB18', peptide_list$proteins) & grepl('ENTR', peptide_list$proteins)), ]

peptide_list$target_L = gsub('I', 'L', peptide_list$targets)
ISB18_target_L = peptide_list$target_L[grepl('ISB18', peptide_list$proteins)]
peptide_list_sub = peptide_list %>% 
  filter(target_L %in% ISB18_target_L) %>% group_by(target_L) %>% mutate(check = if (n() == 1) { FALSE } else { any(grepl('ISB18', proteins)) & any(grepl('ENTR', proteins)) } ) %>% 
  ungroup() %>% filter(check == TRUE) %>% filter(grepl('ENTR', proteins))
peptide_list = peptide_list %>% filter(!(targets %in% peptide_list_sub$targets))

peptide_list$proteins = paste(">", peptide_list$proteins, sep = "") #898 ISB18 peptides, 571298 ENTR peptides

peptide_list = peptide_list[, 2:1]

ISB18_inds = which(grepl('ISB18', peptide_list$proteins))
set.seed(1)
ISB18_inds_remove = sample(ISB18_inds, size = floor(length(ISB18_inds)*3/4))
peptide_list = peptide_list[-ISB18_inds_remove, ] #225 ISB18 peptides, 571298 ENTR peptides

vec_pep_list = as.vector(t(as.matrix(peptide_list)))
file = file("combine_025.fasta")
writeLines(vec_pep_list, file)
close(file)


########################################################################################################################

#creates randomly shuffled databases

#library
library(stringi)
library(seqinr)
options(scipen = 999)

ids = list.files(path = 'datasets')
set.seed(10112023)
for (id in ids) {
  print(id)
  all_files = list.files(path = paste('datasets/', id, sep = ''))
  fasta = all_files[grepl('\\.fasta', all_files)]
  import_fasta = read.fasta(paste('datasets/', id, '/', fasta, sep = ''), forceDNAtolower = FALSE)
  proteins = lapply(import_fasta, function(x) paste(x, collapse = ''))
  
  rev_proteins = lapply(proteins, stri_reverse)
  names_rev_proteins = paste0('rev_', names(rev_proteins))
  names(rev_proteins) = names_rev_proteins
  
  unlink(paste('datasets/', id, '/msfragger_fasta', sep = ''), recursive = TRUE)
  dir.create(file.path(paste('datasets/', id, sep = ''), 'msfragger_fasta'), showWarnings = FALSE)
  
  combined_proteins = c(proteins, rev_proteins)
  write.fasta(combined_proteins, names(combined_proteins), file.out = paste('datasets/', id, '/msfragger_fasta/', 'combined_', 0, '_decoy.fasta', sep = ''))
  
  proteins_split = lapply(proteins, function(x) unlist(strsplit(x, "(?<=[KR])", perl = TRUE)))
  cycle_proteins_split = lapply(proteins_split, function(x) paste(substr(x, nchar(x) - 1, nchar(x) - 1), substr(x, 1, nchar(x) - 2), substr(x, nchar(x), nchar(x)), sep = ""))
  cycle_proteins = lapply(cycle_proteins_split, function(x) paste(x, collapse = ''))
  names_cycle_proteins = paste0('rev_2_', names(cycle_proteins))
  names(cycle_proteins) = names_cycle_proteins
  
  combined_proteins = c(proteins, cycle_proteins)
  write.fasta(combined_proteins, names(combined_proteins), file.out = paste('datasets/', id, '/msfragger_fasta/', 'combined_', 1, '_decoy.fasta', sep = ''))
  
  combined_proteins = c(proteins, rev_proteins, cycle_proteins)
  write.fasta(combined_proteins, names(combined_proteins), file.out = paste('datasets/', id, '/msfragger_fasta/', 'combined', '_decoy.fasta', sep = ''))
}

#warning for PXD012528 but it appears fine on inspection.

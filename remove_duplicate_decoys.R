#runs a script to remove matches that are ambiguous if they are targets or decoys.

datasets = list.files(path = 'datasets')
for (dataset in datasets){
  for (i in 0:1) {
    read_pin = read_delim(paste('datasets/', dataset, '/msfragger_output/narrow_', i, '.pin', sep = ''))
    proteins = str_split(read_pin$Proteins, '\t')
    proteins = lapply(proteins, function(x) x[x!= ''])
    proteins_check = unlist(lapply(proteins, function(x) any(grepl('rev_', x)) & any(!(grepl('rev_', x)))))
    read_pin = read_pin[!proteins_check, ]
    write_delim(read_pin, paste('datasets/', dataset, '/msfragger_output/narrow_dups_removed_', i, '.pin', sep = ''), delim = "\t", quote = 'none')
  }
  
  for (i in 0:1) {
    read_pin = read_delim(paste('datasets/', dataset, '/msfragger_output/open_', i, '.pin', sep = ''))
    proteins = str_split(read_pin$Proteins, '\t')
    proteins = lapply(proteins, function(x) x[x!= ''])
    proteins_check = unlist(lapply(proteins, function(x) any(grepl('rev_', x)) & any(!(grepl('rev_', x)))))
    read_pin = read_pin[!proteins_check, ]
    write_delim(read_pin, paste('datasets/', dataset, '/msfragger_output/open_dups_removed_', i, '.pin', sep = ''), delim = "\t", quote = 'none')
  }
  
  read_pin = read_delim(paste('datasets/', dataset, '/msfragger_output/narrow', '.pin', sep = ''))
  proteins = str_split(read_pin$Proteins, '\t')
  proteins = lapply(proteins, function(x) x[x!= ''])
  proteins_check = unlist(lapply(proteins, function(x) (any(grepl('rev_', x)) & any(!(grepl('rev_', x)))) ))
  read_pin = read_pin[!proteins_check, ]
  write_delim(read_pin, paste('datasets/', dataset, '/msfragger_output/narrow_dups_removed', '.pin', sep = ''), delim = "\t", quote = 'none')
  
  read_pin = read_delim(paste('datasets/', dataset, '/msfragger_output/open', '.pin', sep = ''))
  proteins = str_split(read_pin$Proteins, '\t')
  proteins = lapply(proteins, function(x) x[x!= ''])
  proteins_check = unlist(lapply(proteins, function(x) (any(grepl('rev_', x)) & any(!(grepl('rev_', x)))) ))
  read_pin = read_pin[!proteins_check, ]
  write_delim(read_pin, paste('datasets/', dataset, '/msfragger_output/open_dups_removed', '.pin', sep = ''), delim = "\t", quote = 'none')
}



for (dataset in datasets) {
  for (i in 0:9) {
    peptide_list = read.delim(paste('datasets/', dataset, '/index-', i, '/tide-index.peptides.txt', sep = ''))
    read_pin = read_delim(paste('datasets/', dataset, '/crux-output/narrow_1_', i, '.make-pin.pin', sep = ''))
    remove_aa = unlist(lapply(read_pin$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
    read_pin = read_pin[remove_aa %in% c(peptide_list$target, peptide_list$decoy.s.), ]
    write_delim(read_pin, paste('datasets/', dataset, '/crux-output/narrow_dups_removed_', i, '.make-pin.pin', sep = ''), delim = "\t", quote = 'none')
    
    read_pin = read_delim(paste('datasets/', dataset, '/crux-output/open_1_', i, '.make-pin.pin', sep = ''))
    remove_aa = unlist(lapply(read_pin$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
    read_pin = read_pin[remove_aa %in% c(peptide_list$target, peptide_list$decoy.s.), ]
    write_delim(read_pin, paste('datasets/', dataset, '/crux-output/open_dups_removed_', i, '.make-pin.pin', sep = ''), delim = "\t", quote = 'none')
  }
}




for (dataset in datasets) {
  for (i in 0:9) {
    peptide_list = read.delim(paste('datasets/', dataset, '/index-', i, '/tide-index.peptides.txt', sep = ''))
    read_pin = read_delim(paste('datasets/', dataset, '/crux-output/narrow_5', i, '_sep.tide-search.target.pin', sep = ''))
    remove_aa = unlist(lapply(read_pin$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
    read_pin = read_pin[remove_aa %in% c(peptide_list$target, peptide_list$decoy.s.), ]
    write_delim(read_pin, paste('datasets/', dataset, '/crux-output/narrow_dups_removed_', i, '_sep.tide-search.target.pin', sep = ''), delim = "\t", quote = 'none')
    
    read_pin = read_delim(paste('datasets/', dataset, '/crux-output/narrow_5', i, '_sep.tide-search.decoy.pin', sep = ''))
    remove_aa = unlist(lapply(read_pin$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
    read_pin = read_pin[remove_aa %in% c(peptide_list$target, peptide_list$decoy.s.), ]
    write_delim(read_pin, paste('datasets/', dataset, '/crux-output/narrow_dups_removed_', i, '_sep.tide-search.decoy.pin', sep = ''), delim = "\t", quote = 'none')
    
    
    read_pin = read_delim(paste('datasets/', dataset, '/crux-output/open_5', i, '_sep.tide-search.target.pin', sep = ''))
    remove_aa = unlist(lapply(read_pin$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
    read_pin = read_pin[remove_aa %in% c(peptide_list$target, peptide_list$decoy.s.), ]
    write_delim(read_pin, paste('datasets/', dataset, '/crux-output/open_dups_removed_', i, '_sep.tide-search.target.pin', sep = ''), delim = "\t", quote = 'none')
    
    read_pin = read_delim(paste('datasets/', dataset, '/crux-output/open_5', i, '_sep.tide-search.decoy.pin', sep = ''))
    remove_aa = unlist(lapply(read_pin$Peptide, function(x) substr(x, 3, nchar(x) - 2)))
    read_pin = read_pin[remove_aa %in% c(peptide_list$target, peptide_list$decoy.s.), ]
    write_delim(read_pin, paste('datasets/', dataset, '/crux-output/open_dups_removed_', i, '_sep.tide-search.decoy.pin', sep = ''), delim = "\t", quote = 'none')
    
  }
}


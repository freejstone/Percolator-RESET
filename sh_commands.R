command = "crux tide-index --overwrite T --seed IND --auto-modifications T --auto-modifications-spectra 'SPECTRUM_FILE' --peptide-list T --output-dir index-IND 'FASTA_FILE' index-IND"
datasets = list.files(path = 'datasets', pattern = 'PXD')
command_all = c()
for (dataset in datasets){
  files = list.files(path = paste('datasets/', dataset, sep = ''))
  spectrum_file = files[grepl('mgf', files, ignore.case = TRUE)]
  fasta_file = files[grepl('fasta', files, ignore.case = TRUE)]
  command_all = c(command_all, paste('cd ', dataset, sep = ''))
  for (i in 0:9) {
    command_temp = gsub('IND', i, command)
    command_temp = gsub('SPECTRUM_FILE', spectrum_file, command_temp)
    command_temp = gsub('FASTA_FILE', fasta_file, command_temp)
    command_all = c(command_all, command_temp)
  }
  command_all = c(command_all, 'cd ..')
}

command = "crux tide-search --auto-precursor-window warn --auto-mz-bin-width warn --use-tailor-calibration T --concat T --top-match 5 --num-threads 30 --overwrite T --fileroot narrow_5IND 'SPECTRUM_FILE' index-IND"
datasets = list.files(path = 'datasets', pattern = 'PXD')
command_all = c()
for (dataset in datasets){
  files = list.files(path = paste('datasets/', dataset, sep = ''))
  spectrum_file = files[grepl('mgf', files, ignore.case = TRUE)]
  command_all = c(command_all, paste('cd ', dataset, sep = ''))
  for (i in 0:9) {
    command_temp = gsub('IND', i, command)
    command_temp = gsub('SPECTRUM_FILE', spectrum_file, command_temp)
    command_all = c(command_all, command_temp)
  }
  command_all = c(command_all, 'cd ..')
}

command = "crux tide-search --auto-mz-bin-width warn --precursor-window-type mass --precursor-window 100 --use-tailor-calibration T --concat T --top-match 5 --num-threads 4 --overwrite T --fileroot open_5IND 'SPECTRUM_FILE' index-IND"
datasets = list.files(path = 'datasets', pattern = 'PXD')
command_all = c()
for (dataset in datasets){
  files = list.files(path = paste('datasets/', dataset, sep = ''))
  spectrum_file = files[grepl('mgf', files, ignore.case = TRUE)]
  command_all = c(command_all, paste('cd ', dataset, sep = ''))
  for (i in 0:9) {
    command_temp = gsub('IND', i, command)
    command_temp = gsub('SPECTRUM_FILE', spectrum_file, command_temp)
    command_all = c(command_all, command_temp)
  }
  command_all = c(command_all, 'cd ..')
}


command = "crux make-pin --overwrite T --fileroot search_file_root crux-output/search_file"
datasets = list.files(path = 'datasets', pattern = 'PXD')
command_all = c()
for (dataset in datasets){
  command_all = c(command_all, paste('cd ', dataset, sep = ''))
  for (i in 0:9){
    command_temp = gsub('search_file_root', paste('open_5', i, sep = ''), command)
    command_temp = gsub('search_file', paste('open_5', i, '.tide-search.txt', sep = ''), command_temp)
    command_all = c(command_all, command_temp)
    
    command_temp = gsub('search_file_root', paste('narrow_5', i, sep = ''), command)
    command_temp = gsub('search_file', paste('narrow_5', i, '.tide-search.txt', sep = ''), command_temp)
    command_all = c(command_all, command_temp)
  }
  command_all = c(command_all, 'cd ..')
}



command = "python do_super_percolator.py --seed 0 --overwrite T --output_dir datasets/PXID/crux-output --file_root IND datasets/PXID/crux-output/narrow_5IND.make-pin.pin datasets/PXID/crux-output/open_5IND.make-pin.pin datasets/PXID/index-IND/tide-index.peptides.txt"
datasets = list.files(path = 'datasets', pattern = 'PXD')
command_all = c()
for (dataset in datasets){
  for (i in 0:9){
    command_temp = gsub('PXID', dataset, command)
    command_temp = gsub('IND', i, command_temp)
    command_all = c(command_all, command_temp)
  }
}

command = "crux percolator --overwrite T --tdc F --only-psms T --output-dir datasets/PXID/crux-output --fileroot open_1_IND datasets/PXID/crux-output/open_1_IND.make-pin.pin"
datasets = list.files(path = 'datasets', pattern = 'PXD')
command_all = c()
for (dataset in datasets){
  for (i in 0:9){
    command_temp = gsub('PXID', dataset, command)
    command_temp = gsub('IND', i, command_temp)
    command_all = c(command_all, command_temp)
  }
}


command = "python do_super_percolator.py --seed 0 --overwrite T --top_positive T --train_FDR_threshold 0.005 --output_dir datasets/PXID/crux-output --file_root IND_top_positive datasets/PXID/crux-output/narrow_5IND.make-pin.pin datasets/PXID/crux-output/open_5IND.make-pin.pin datasets/PXID/index-IND/tide-index.peptides.txt"
datasets = list.files(path = 'datasets', pattern = 'PXD')
command_all = c()
for (dataset in datasets){
  for (i in 0:9){
    command_temp = gsub('PXID', dataset, command)
    command_temp = gsub('IND', i, command_temp)
    command_all = c(command_all, command_temp)
  }
}

command = "python do_super_percolator.py --seed 0 --overwrite T --top_positive T --train_FDR_threshold 0.05 --svm F --output_dir datasets/PXID/crux-output --file_root IND_qda datasets/PXID/crux-output/narrow_5IND.make-pin.pin datasets/PXID/crux-output/open_5IND.make-pin.pin datasets/PXID/index-IND/tide-index.peptides.txt"
datasets = list.files(path = 'datasets', pattern = 'PXD')
command_all = c()
for (dataset in datasets){
  for (i in 0:9){
    command_temp = gsub('PXID', dataset, command)
    command_temp = gsub('IND', i, command_temp)
    command_all = c(command_all, command_temp)
  }
}



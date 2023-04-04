#Compare to percolator

setwd("/Volumes/Seagate Backup Plus Drive/JACK/PhD/ppx_files/High resolution/Static_modification/PXD016724/crux-output/super_percolator/search-0")

#functions
TDC_flex_c = function(decoy_wins, target_wins, BC1 = 1, c = 1/2, lambda = 1/2){
  #give in order of best score first to worst score last
  nTD <- cumsum(target_wins)
  nDD <- cumsum(decoy_wins)
  fdps <- (pmin(1, ((BC1 + nDD)/ nTD) * (c / (1-lambda)))) #nDD, nTD have same length
  qvals <- rev(cummin(rev(fdps)))
  return(qvals)
}

#files
target = read.delim('./crux-output/ord.percolator.target.psms.txt') #19276
decoy = read.delim('./crux-output/ord.percolator.decoy.psms.txt') #11841

#do TDC on percolator scores
target_decoy = rbind(target, decoy)
target_decoy = target_decoy[sample(nrow(target_decoy)), ]
target_decoy = target_decoy[order(-target_decoy$score), ]
target_decoy$peptide = gsub(".*[.]([^.]+)[.].*", "\\1", target_decoy$peptide)
target_decoy = target_decoy[!duplicated(target_decoy$peptide), ]

q_val = TDC_flex_c(grepl('decoy', target_decoy$PSMId), grepl('target', target_decoy$PSMId))
sum(q_val <= 0.01 & grepl('target', target_decoy$PSMId)) #gives 4736 discoveries

#files
target = read.delim('./crux-output/percolator.target.psms.txt') #31108
decoy = read.delim('./crux-output/percolator.decoy.psms.txt') #31108

target_decoy = rbind(target, decoy)
target_decoy = target_decoy[sample(nrow(target_decoy)), ]
target_decoy = target_decoy[order(-target_decoy$score), ]
target_decoy$peptide = gsub(".*[.]([^.]+)[.].*", "\\1", target_decoy$peptide)
target_decoy = target_decoy[!duplicated(target_decoy$peptide), ]

q_val = TDC_flex_c(grepl('decoy', target_decoy$PSMId), grepl('target', target_decoy$PSMId))
sum(q_val <= 0.01 & grepl('target', target_decoy$PSMId)) #gives 4171 discoveries

target$target_decoy = 'target'
decoy$target_decoy = 'decoy'
target$PSMId = gsub(".*target[_]0[_]([^.]+)[_]+.*[_][^.]+.*", "\\1", target$PSMId)
decoy$PSMId = gsub(".*decoy[_]0[_]([^.]+)[_]+.*[_][^.]+.*", "\\1", decoy$PSMId)

target_decoy = rbind(target, decoy)
target_decoy = target_decoy[sample(nrow(target_decoy)), ]
target_decoy = target_decoy[order(-target_decoy$score), ]
target_decoy = target_decoy[!duplicated(target_decoy$PSMId), ]
sum(target_decoy$target_decoy == 'target') #19263
sum(target_decoy$target_decoy == 'decoy') #11845
target_decoy$peptide = gsub(".*[.]([^.]+)[.].*", "\\1", target_decoy$peptide)
target_decoy = target_decoy[!duplicated(target_decoy$peptide), ]

q_val = TDC_flex_c(grepl('decoy', target_decoy$target_decoy), grepl('target', target_decoy$target_decoy))
sum(q_val <= 0.01 & grepl('target', target_decoy$target_decoy)) #gives 4709 discoveries

setwd("/Volumes/Seagate Backup Plus Drive/JACK/PhD/ppx_files/High resolution/Static_modification/PXD016724/crux-output/super_percolator/search-concat-0")

#files
target = read.delim('./crux-output/ord.percolator.target.psms.txt') #19270
decoy = read.delim('./crux-output/ord.percolator.decoy.psms.txt') #11838

peptide_list_0 = read.delim('../index-0/tide-index.peptides.txt')

target$target_decoy = 'target'
decoy$target_decoy = 'decoy'
target$peptide = gsub(".*[.]([^.]+)[.].*", "\\1", target$peptide)
decoy$peptide = gsub(".*[.]([^.]+)[.].*", "\\1", decoy$peptide)

target$original = target$peptide
indxs = match(decoy$peptide, peptide_list_0$decoy)
decoy$original = peptide_list_0$target[indxs]

target_decoy = rbind(target, decoy)
target_decoy = target_decoy[sample(nrow(target_decoy)), ]
target_decoy = target_decoy[order(-target_decoy$score), ]
target_decoy = target_decoy[!duplicated(target_decoy$original), ]

q_val = TDC_flex_c(grepl('decoy', target_decoy$target_decoy), grepl('target', target_decoy$target_decoy))
sum(q_val <= 0.01 & grepl('target', target_decoy$target_decoy)) #4821 (peptide-level competition is made redundant somewhat)


library(dplyr)

TDC_flex_c = function(decoy_wins, target_wins, BC1 = 1, c = 1/2, lambda = 1/2){
  #give in order of best score first to worst score last
  nTD = cumsum(target_wins)
  nDD = cumsum(decoy_wins)
  fdps = (pmin(1, ((BC1 + nDD)/ nTD) * (c / (1-lambda)))) #nDD, nTD have same length
  qvals = rev(cummin(rev(fdps)))
  return(qvals)
}
do_TDC = function(target_decoys_TDC){
  target_decoys_TDC = target_decoys_TDC[sample(nrow(target_decoys_TDC)), ]
  #take the best scoring PSM according to cluster + original.target.sequence
  target_decoys_TDC = target_decoys_TDC[order(-target_decoys_TDC$score), ]
  
  target_decoys_TDC = target_decoys_TDC[ !duplicated(target_decoys_TDC[, c('original_target_sequence')]), ]
  
  #get winning information
  winning_scores = target_decoys_TDC$score
  
  labels = target_decoys_TDC$target.decoy
  labels[labels == 'target'] = 1
  labels[labels == 'decoy'] = -1
  labels = as.numeric(labels)
  
  winning_peptides = target_decoys_TDC$peptide
  
  df = data.frame(winning_scores, winning_peptides, labels)
  
  df = df[sample(1:nrow(df)),]
  
  df = df[order(-df$winning_scores), ]
  return(df)
}

df_all = data.frame(matrix(0, ncol = 4))
colnames(df_all) = c("Threshold", "power_open", "dcy_indx", "PXID")
alpha_vec = seq(0.01, 0.05, by = 0.001)
d_vec = c(0:9)
count_all = 1
set.seed(1)
datasets = list.files(path = 'datasets', pattern = 'PXD')

for (dataset in datasets) {
  print(dataset)
  for (d in d_vec) {
    txt = paste('datasets', '/', dataset, '/index-', d, '/tide-index.peptides.txt', sep = "")
    peptide_list = read.delim(txt)
    
    txt_decoy = paste('datasets', '/', dataset, '/crux-output/open_1_', d, '.percolator.decoy.psms.txt', sep = "")
    txt_target = paste('datasets', '/', dataset, '/crux-output/open_1_', d, '.percolator.target.psms.txt', sep = "")
    
    targets = read.delim(txt_target)
    targets$target.decoy = 'target'
    
    decoys = read.delim(txt_decoy)
    decoys$target.decoy = 'decoy'
    
    targets$original_target_sequence = substring(targets$peptide, 3, nchar(targets$peptide) - 2)
    targets$original_target_sequence = gsub('[0-9]+|\\.|\\[|\\]', '', targets$original_target_sequence)
    
    decoys$original_target_sequence = substring(decoys$peptide, 3, nchar(decoys$peptide) - 2)
    decoys$original_target_sequence = gsub('[0-9]+|\\.|\\[|\\]', '', decoys$original_target_sequence)
    decoys$original_target_sequence = peptide_list$target[match(decoys$original_target_sequence, peptide_list$decoy)]
    
    target_decoys = rbind(targets, decoys)
    
    df_open = do_TDC(target_decoys)
    
    q_vals_open = TDC_flex_c(df_open$labels == -1, df_open$labels == 1)
    
    ##############################################################################################################################
    
    for (alpha in alpha_vec){
      df_all[count_all, ] <- c(alpha, sum(q_vals_open <= alpha & df_open$labels == 1), d, dataset)
      count_all <- count_all + 1
    }
  }
}


write.csv(df_all, paste('percolator_all.csv', sep = ''))

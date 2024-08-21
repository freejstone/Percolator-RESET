#functions
TDC_flex_c = function(decoy_wins, target_wins, BC1 = 1, c = 1/2, lambda = 1/2){
  #give in order of best score first to worst score last
  nTD = cumsum(target_wins)
  nDD = cumsum(decoy_wins)
  fdps = (pmin(1, ((BC1 + nDD)/ nTD) * (c / (1-lambda)))) #nDD, nTD have same length
  qvals = rev(cummin(rev(fdps)))
  return(qvals)
}

#library
library(tidyverse)

df_all = data.frame(matrix(0, ncol = 5))
colnames(df_all) = c("Threshold", "power", "est_FDP", "index", "sample_sizes")
alpha_vec =seq(0.01, 0.1, by = 0.001)
d_vec = c(0:99)
sample_sizes = c('full', '075', 'half', '025')
seed = 1
count_all = 1
set.seed(1234)
for (s in sample_sizes) {
  for (d in d_vec){
    
    print(d)
    txt = paste('percolator_results/', 'open_top1_enzint_zero_', s, '_', d,'.percolator.target.psms.txt', sep = "")
    targets = read.delim(txt)
    txt = paste('percolator_results/', 'open_top1_enzint_zero_', s, '_', d,'.percolator.decoy.psms.txt', sep = "")
    decoys = read.delim(txt)
    
    target_decoys = rbind(targets, decoys)
    
    txt = paste('./tide-index/index-', s, '-', d, '/tide-index.peptides.txt', sep = "")
    
    peptide_list = read.delim(txt)
    names(peptide_list)[2] = 'decoy'
    
    target_decoys = target_decoys[sample(nrow(target_decoys)), ]
    
    target_decoys$target_decoy = grepl('decoy_', target_decoys$PSMId)
    target_decoys$target_decoy[target_decoys$target_decoy == 'TRUE'] = 'decoy'
    target_decoys$target_decoy[target_decoys$target_decoy == 'FALSE'] = 'target'
    
    target_decoys$peptide = unlist(lapply(target_decoys$peptide, function(x) unlist(strsplit(x, "\\."))[2]))
    
    check = all(target_decoys$peptide %in% c(peptide_list$target, peptide_list$decoy))
    print(check)
    if (!check){
      print('OH NOOO')
    }
    
    target_decoys$original_target_sequence = target_decoys$peptide
    target_decoys$original_target_sequence[target_decoys$target_decoy == 'decoy'] = peptide_list$target[match(target_decoys$peptide[target_decoys$target_decoy == 'decoy'], peptide_list$decoy)]
    
    target_decoys = target_decoys[order(-target_decoys$score), ]
    
    target_decoys = target_decoys[!duplicated(target_decoys$original_target_sequence), ]
    
    target_decoys$proteinIds = peptide_list$proteins[match(target_decoys$original_target_sequence, peptide_list$target)]
    
    winning_score = target_decoys$score
    winning_protein = target_decoys$proteinIds
    labels = target_decoys$target_decoy
    labels[labels == 'decoy'] = -1
    labels[labels == 'target'] = 1
    
    df = data.frame(winning_score, winning_protein, labels)
    
    df = df[sample(nrow(df)), ]
    
    df = df[order(-df$winning_score), ]
    
    q_vals = TDC_flex_c(df$labels == -1, df$labels == 1)
    for (alpha in alpha_vec){
      power = sum(q_vals <= alpha & df$labels == 1)
      proteins_identified = df$winning_protein[q_vals <= alpha & df$labels == 1]
      
      est_FDP = sum(grepl('ENTR', proteins_identified))/max(power, 1)
      if (alpha == 0.01){
        print(est_FDP)
      }
      df_all[count_all, ] = c(alpha, power, est_FDP, d, s)
      count_all = count_all + 1
    }
  }
}

write.csv(df_all, 'results/percolator_open_tide.csv')


df_all = data.frame(matrix(0, ncol = 5))
colnames(df_all) = c("Threshold", "power", "est_FDP", "index", "sample_sizes")
alpha_vec =seq(0.01, 0.1, by = 0.001)
d_vec = c(0:99)
sample_sizes = c('full', '075', 'half', '025')
seed = 1
count_all = 1
set.seed(1234)
for (s in sample_sizes) {
  for (d in d_vec){
    
    print(d)
    txt = paste('percolator_results/', 'dcy_enzint_zero_', s, '_', d,'.percolator.target.psms.txt', sep = "")
    targets = read.delim(txt)
    txt = paste('percolator_results/', 'dcy_enzint_zero_', s, '_', d,'.percolator.decoy.psms.txt', sep = "")
    decoys = read.delim(txt)
    
    target_decoys = rbind(targets, decoys)
    
    txt = paste('./tide-index/index-', s, '-', d, '/tide-index.peptides.txt', sep = "")
    
    peptide_list = read.delim(txt)
    names(peptide_list)[2] = 'decoy'
    
    target_decoys = target_decoys[sample(nrow(target_decoys)), ]
    
    target_decoys$target_decoy = grepl('decoy_', target_decoys$PSMId)
    target_decoys$target_decoy[target_decoys$target_decoy == 'TRUE'] = 'decoy'
    target_decoys$target_decoy[target_decoys$target_decoy == 'FALSE'] = 'target'
    
    target_decoys$peptide = unlist(lapply(target_decoys$peptide, function(x) unlist(strsplit(x, "\\."))[2]))
    
    check = all(target_decoys$peptide %in% c(peptide_list$target, peptide_list$decoy))
    print(check)
    if (!check){
      print('OH NOOO')
    }
    
    target_decoys$original_target_sequence = target_decoys$peptide
    target_decoys$original_target_sequence[target_decoys$target_decoy == 'decoy'] = peptide_list$target[match(target_decoys$peptide[target_decoys$target_decoy == 'decoy'], peptide_list$decoy)]
    
    target_decoys = target_decoys[order(-target_decoys$score), ]
    
    target_decoys = target_decoys[!duplicated(target_decoys$original_target_sequence), ]
    
    target_decoys$proteinIds = peptide_list$proteins[match(target_decoys$original_target_sequence, peptide_list$target)]
    
    winning_score = target_decoys$score
    winning_protein = target_decoys$proteinIds
    labels = target_decoys$target_decoy
    labels[labels == 'decoy'] = -1
    labels[labels == 'target'] = 1
    
    df = data.frame(winning_score, winning_protein, labels)
    
    df = df[sample(nrow(df)), ]
    
    df = df[order(-df$winning_score), ]
    
    q_vals = TDC_flex_c(df$labels == -1, df$labels == 1)
    for (alpha in alpha_vec){
      power = sum(q_vals <= alpha & df$labels == 1)
      proteins_identified = df$winning_protein[q_vals <= alpha & df$labels == 1]
      
      est_FDP = sum(grepl('ENTR', proteins_identified))/max(power, 1)
      if (alpha == 0.01){
        print(est_FDP)
      }
      df_all[count_all, ] = c(alpha, power, est_FDP, d, s)
      count_all = count_all + 1
    }
  }
}

write.csv(df_all, 'results/percolator_narrow_tide.csv')


############################################################



#functions
TDC_flex_c = function(decoy_wins, target_wins, BC1 = 1, c = 1/2, lambda = 1/2){
  #give in order of best score first to worst score last
  nTD = cumsum(target_wins)
  nDD = cumsum(decoy_wins)
  fdps = (pmin(1, ((BC1 + nDD)/ nTD) * (c / (1-lambda)))) #nDD, nTD have same length
  qvals = rev(cummin(rev(fdps)))
  return(qvals)
}

#library
library(tidyverse)

df_all = data.frame(matrix(0, ncol = 5))
colnames(df_all) = c("Threshold", "power", "est_FDP", "index", "sample_sizes")
alpha_vec =seq(0.01, 0.1, by = 0.001)
d_vec = c(0:99)
sample_sizes = c('full', '075', 'half', '025')
seed = 1
count_all = 1
set.seed(1234)
for (s in sample_sizes) {
  for (d in d_vec){
    
    print(d)
    txt = paste('percolator_results/', 'open_top1_enzint_zero_', s, '_', d,'.percolator.target.psms.txt', sep = "")
    targets = read.delim(txt)
    txt = paste('percolator_results/', 'open_top1_enzint_zero_', s, '_', d,'.percolator.decoy.psms.txt', sep = "")
    decoys = read.delim(txt)
    
    target_decoys = rbind(targets, decoys)
    
    txt = paste('./tide-index/index-', s, '-', d, '/tide-index.peptides.txt', sep = "")
    
    peptide_list = read.delim(txt)
    names(peptide_list)[2] = 'decoy'
    
    target_decoys = target_decoys[sample(nrow(target_decoys)), ]
    
    target_decoys$target_decoy = grepl('decoy_', target_decoys$PSMId)
    target_decoys$target_decoy[target_decoys$target_decoy == 'TRUE'] = 'decoy'
    target_decoys$target_decoy[target_decoys$target_decoy == 'FALSE'] = 'target'
    
    target_decoys$peptide = unlist(lapply(target_decoys$peptide, function(x) unlist(strsplit(x, "\\."))[2]))
    
    check = all(target_decoys$peptide %in% c(peptide_list$target, peptide_list$decoy))
    print(check)
    if (!check){
      print('OH NOOO')
    }
    
    target_decoys$original_target_sequence = target_decoys$peptide
    target_decoys$original_target_sequence[target_decoys$target_decoy == 'decoy'] = peptide_list$target[match(target_decoys$peptide[target_decoys$target_decoy == 'decoy'], peptide_list$decoy)]
    
    target_decoys = target_decoys[order(-target_decoys$score), ]
    
    target_decoys = target_decoys[!duplicated(target_decoys$peptide), ]
    
    target_decoys$proteinIds = peptide_list$proteins[match(target_decoys$original_target_sequence, peptide_list$target)]
    
    winning_score = target_decoys$score
    winning_protein = target_decoys$proteinIds
    labels = target_decoys$target_decoy
    labels[labels == 'decoy'] = -1
    labels[labels == 'target'] = 1
    
    df = data.frame(winning_score, winning_protein, labels)
    
    df = df[sample(nrow(df)), ]
    
    df = df[order(-df$winning_score), ]
    
    q_vals = TDC_flex_c(df$labels == -1, df$labels == 1)
    for (alpha in alpha_vec){
      power = sum(q_vals <= alpha & df$labels == 1)
      proteins_identified = df$winning_protein[q_vals <= alpha & df$labels == 1]
      
      est_FDP = sum(grepl('ENTR', proteins_identified))/max(power, 1)
      if (alpha == 0.01){
        print(est_FDP)
      }
      df_all[count_all, ] = c(alpha, power, est_FDP, d, s)
      count_all = count_all + 1
    }
  }
}

write.csv(df_all, 'results/percolator_open_tide_psm_only.csv')


df_all = data.frame(matrix(0, ncol = 5))
colnames(df_all) = c("Threshold", "power", "est_FDP", "index", "sample_sizes")
alpha_vec =seq(0.01, 0.1, by = 0.001)
d_vec = c(0:99)
sample_sizes = c('full', '075', 'half', '025')
seed = 1
count_all = 1
set.seed(1234)
for (s in sample_sizes) {
  for (d in d_vec){
    
    print(d)
    txt = paste('percolator_results/', 'dcy_enzint_zero_', s, '_', d,'.percolator.target.psms.txt', sep = "")
    targets = read.delim(txt)
    txt = paste('percolator_results/', 'dcy_enzint_zero_', s, '_', d,'.percolator.decoy.psms.txt', sep = "")
    decoys = read.delim(txt)
    
    target_decoys = rbind(targets, decoys)
    
    txt = paste('./tide-index/index-', s, '-', d, '/tide-index.peptides.txt', sep = "")
    
    peptide_list = read.delim(txt)
    names(peptide_list)[2] = 'decoy'
    
    target_decoys = target_decoys[sample(nrow(target_decoys)), ]
    
    target_decoys$target_decoy = grepl('decoy_', target_decoys$PSMId)
    target_decoys$target_decoy[target_decoys$target_decoy == 'TRUE'] = 'decoy'
    target_decoys$target_decoy[target_decoys$target_decoy == 'FALSE'] = 'target'
    
    target_decoys$peptide = unlist(lapply(target_decoys$peptide, function(x) unlist(strsplit(x, "\\."))[2]))
    
    check = all(target_decoys$peptide %in% c(peptide_list$target, peptide_list$decoy))
    print(check)
    if (!check){
      print('OH NOOO')
    }
    
    target_decoys$original_target_sequence = target_decoys$peptide
    target_decoys$original_target_sequence[target_decoys$target_decoy == 'decoy'] = peptide_list$target[match(target_decoys$peptide[target_decoys$target_decoy == 'decoy'], peptide_list$decoy)]
    
    target_decoys = target_decoys[order(-target_decoys$score), ]
    
    target_decoys = target_decoys[!duplicated(target_decoys$peptide), ]
    
    target_decoys$proteinIds = peptide_list$proteins[match(target_decoys$original_target_sequence, peptide_list$target)]
    
    winning_score = target_decoys$score
    winning_protein = target_decoys$proteinIds
    labels = target_decoys$target_decoy
    labels[labels == 'decoy'] = -1
    labels[labels == 'target'] = 1
    
    df = data.frame(winning_score, winning_protein, labels)
    
    df = df[sample(nrow(df)), ]
    
    df = df[order(-df$winning_score), ]
    
    q_vals = TDC_flex_c(df$labels == -1, df$labels == 1)
    for (alpha in alpha_vec){
      power = sum(q_vals <= alpha & df$labels == 1)
      proteins_identified = df$winning_protein[q_vals <= alpha & df$labels == 1]
      
      est_FDP = sum(grepl('ENTR', proteins_identified))/max(power, 1)
      if (alpha == 0.01){
        print(est_FDP)
      }
      df_all[count_all, ] = c(alpha, power, est_FDP, d, s)
      count_all = count_all + 1
    }
  }
}

write.csv(df_all, 'results/percolator_narrow_tide_psm_only.csv')


############################################################


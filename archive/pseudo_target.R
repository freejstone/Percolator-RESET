#setting wd
setwd("~/Documents/PhD/ppx_files/High resolution/Static_modification/PXD025130/crux-output/super_percolator")

#setting options and installing packages
options(scipen = 999)
library(tidyverse)
library(caret)
library(e1071)

#functions
TDC_flex_c = function(decoy_wins, target_wins, BC1 = 1, c = 1/2, lambda = 1/2){
  #give in order of best score first to worst score last
  nTD <- cumsum(target_wins)
  nDD <- cumsum(decoy_wins)
  fdps <- (pmin(1, ((BC1 + nDD)/ nTD) * (c / (1-lambda)))) #nDD, nTD have same length
  qvals <- rev(cummin(rev(fdps)))
  return(qvals)
}
radial_kernel = function(x, y) {
  dim = length(x) + 1
  return(exp(-rowSums((x - y)^2)/dim))
}

#taking enzint zero
for (search in c('search-0', 'search-1')) {
  for (td in c('target', 'decoy')) {
    search_output = read_delim(paste('./', search, '/tide-search.', td, '.pin', sep = ''), delim = '\t')
    search_output$enzInt = 0
    write.table(search_output, paste('./', search, '/enzint_zero.tide-search.', td, '.pin', sep = ''), row.names = F, sep = '\t', quote = F)
  }
}

#read in percolator file from search-0
target_file = read_delim('./search-0/enzint_zero.tide-search.target.pin', delim = '\t')
decoy_file = read_delim('./search-0/enzint_zero.tide-search.decoy.pin', delim = '\t')

#drop aa before and after peptide
target_file$Peptide = gsub(".*[.]([^.]+)[.].*", "\\1", target_file$Peptide)
decoy_file$Peptide = gsub(".*[.]([^.]+)[.].*", "\\1", decoy_file$Peptide)

#spectrum level competition
decoy_file = decoy_file[match(target_file$ScanNr, decoy_file$ScanNr), ] 
greater_than = which(target_file$TailorScore > decoy_file$TailorScore)
ties = which(target_file$TailorScore == decoy_file$TailorScore)
greater_than_inds = sample(c(0, 1), size = length(ties), replace = T)
greater_than = c(greater_than, ties[as.logical(greater_than_inds)])
combined_file = decoy_file
combined_file[greater_than, ] = target_file[greater_than, ]

#getting best match to each peptide
combined_file = combined_file[order(-combined_file$TailorScore), ]
combined_file = combined_file[!duplicated(combined_file$Peptide), ]

#reading second decoy file
decoy_second_file = read_delim('./search-1/enzint_zero.tide-search.decoy.pin', delim = '\t')
decoy_second_file$Peptide = gsub(".*[.]([^.]+)[.].*", "\\1", decoy_second_file$Peptide)

#getting best match to each peptide
decoy_second_file = decoy_second_file[order(-decoy_second_file$TailorScore), ]
decoy_second_file = decoy_second_file[!duplicated(decoy_second_file$Peptide), ]

#reading in peptide-list
peptide_list_0 = read_delim('./index-0/tide-index.peptides.txt')
peptide_list_1 = read_delim('./index-1/tide-index.peptides.txt')

if (all(peptide_list_0$target == peptide_list_1$target)) {
  peptide_list = data.frame(peptide_list_0$target, peptide_list_0$`decoy(s)`, peptide_list_1$`decoy(s)`)
  names(peptide_list) = c('target', 'decoy1', 'decoy2')
  
  #get original targets
  combined_file$target = combined_file$Peptide
  indxs = match(combined_file$Peptide[combined_file$Label == -1], peptide_list$decoy1)
  combined_file$target[combined_file$Label == -1] = peptide_list$target[indxs]
  
  #peptide-level competition
  combined_file = combined_file[order(-combined_file$TailorScore), ]
  combined_file = combined_file[!duplicated(combined_file$target), ]
  
  #get original targets for second decoy file
  indxs = match(decoy_second_file$Peptide, peptide_list$decoy2)
  decoy_second_file$target = peptide_list$target[indxs]
  
  #reduce peptide list such that it only contains detected targets
  peptide_list = peptide_list[peptide_list$target %in% combined_file$target, ]
  
  #reduce decoy_second_file so only decoys of detected targets
  #this removes some decoys, maybe we should keep them?
  decoy_second_file = decoy_second_file[decoy_second_file$target %in% combined_file$target, ]
  
  #reorder decoy_second_file to match combined_file
  #not that many pairs, might need to go to more than top 1 rank
  combine_file_sub = combined_file[combined_file$target %in% decoy_second_file$target, ]
  indxs = match(combine_file_sub$target, decoy_second_file$target)
  decoy_second_file = decoy_second_file[indxs, ]
  
  pseudo_winning = decoy_second_file
  winning_inds = which(combine_file_sub$TailorScore > decoy_second_file$TailorScore)
  ties = which(combine_file_sub$TailorScore == decoy_second_file$TailorScore)
  winning_inds = c(winning_inds, ties[as.logical(sample(c(0, 1), length(ties), replace = TRUE))])
  pseudo_winning[winning_inds, ] = combine_file_sub[winning_inds, ]
  pseudo_winning$pseudo_Label = -1
  pseudo_winning$pseudo_Label[winning_inds] = 1
  pseudo_winning = pseudo_winning[order(-pseudo_winning$TailorScore), ]
  
  pseudo_winning$pseudo_Label = factor(pseudo_winning$pseudo_Label)
  pseudo_winning_svm = pseudo_winning %>% select(-c(Peptide, Proteins, target, SpecId)) #remove these for training
  pseudo_winning_svm = pseudo_winning_svm %>% select(-c(Label, ScanNr, deltLCn, Charge1, enzN, enzC, enzInt)) #these have zero variance
  #note that we've expressed charge as a dummy vector via pin file which we can treat now as a numerical
  
  pseudo_winning_svm_scale = pseudo_winning_svm %>% select(-c(pseudo_Label)) %>% scale()
  pseudo_winning_svm_scale = data.frame(pseudo_winning_svm_scale, pseudo_winning_svm$pseudo_Label)
  pseudo_winning_svm_scale = pseudo_winning_svm_scale[order(-pseudo_winning_svm_scale$TailorScore), ]
  
  qvals = TDC_flex_c(pseudo_winning_svm_scale$pseudo_winning_svm.pseudo_Label == -1, pseudo_winning_svm_scale$pseudo_winning_svm.pseudo_Label == 1)
  pseudo_winning_svm_scale$svm_labels = -1
  pseudo_winning_svm_scale$svm_labels[qvals <= 0.05 & pseudo_winning_svm_scale$pseudo_winning_svm.pseudo_Label == 1] = 1
  svm_labels = pseudo_winning_svm_scale$svm_labels
  pseudo_winning_svm_scale = pseudo_winning_svm_scale %>% select(-c(pseudo_winning_svm.pseudo_Label, svm_labels))
  
  costs = c(10, 50, 100, 500, 1000, 5000, 10e3, 50e3, 10e4, 50e4, 10e5)
  power_ests = rep(0, length(costs))
  power_max = 0
  # 
  # count = 1
  # for (cost in costs) {
  #   classifier = svm(formula = pseudo_winning_svm$pseudo_Label ~ .,
  #                    data = pseudo_winning_svm_scale,
  #                    type = 'C-classification',
  #                    kernel = 'radial',
  #                    cost = cost,
  #                    tolerance = 0.1)
  #   pseudo_winning_svm_train = pseudo_winning_svm
  #   pseudo_winning_svm_train$svm_score = classifier$decision.values
  #   pseudo_winning_svm_train = pseudo_winning_svm_train[order(-pseudo_winning_svm_train$svm_score), ]
  #   qval = TDC_flex_c(pseudo_winning_svm_train$pseudo_Label == -1, pseudo_winning_svm_train$pseudo_Label == 1)
  #   power = sum(qval <= 0.01 & pseudo_winning_svm_train$pseudo_Label == 1)
  #   if (power >= power_max) {
  #     classifier_max = classifier
  #   }
  #   power_ests[count] = power
  #   count = count + 1
  #   print(c(cost, power))
  # }
  # 
  # combined_file_features = combined_file %>% select(-c(Peptide, Proteins, target, SpecId))
  # combined_file_features = combined_file_features %>% select(-c(Label, ScanNr, deltLCn, Charge1, enzN, enzC, enzInt))
  # combined_file_features_scale = combined_file_features %>% scale()
  # combined_file_features = data.frame(combined_file_features_scale)
  # 
  # predicted_vals = t(classifier_max$coefs)%*%t(apply(classifier_max$SV, 1, radial_kernel, y = combined_file_features)) - classifier$rho
  # combined_file$svm_score = as.vector(predicted_vals)
  # combined_file = combined_file[order(-combined_file$svm_score), ]
  # qval = TDC_flex_c(combined_file$Label == -1, combined_file$Label == 1)
  # sum(qval <= 0.01 & combined_file$Label == 1)
}


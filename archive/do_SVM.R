#libraries
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
do_iterative_svm = function(target_file_name_1, decoy_file_name_1, decoy_file_name_2, peptide_list_name, FDR_threshold = 0.01, Cs = c(0.1, 1, 10), total_iter = 10, kernel = 'linear') {
  
  target_file_1 = read_delim(target_file_name_1, na = character())
  decoy_file_1 = read_delim(decoy_file_name_1, na = character())
  decoy_file_2 = read_delim(decoy_file_name_2, na = character())
  peptide_list = read_delim(peptide_list_name, na = character())
  
  #fix peptide format
  target_file_1$Peptide = gsub(".*[.]([^.]+)[.].*", "\\1", target_file_1$Peptide)
  decoy_file_1$Peptide = gsub(".*[.]([^.]+)[.].*", "\\1", decoy_file_1$Peptide)
  decoy_file_2$Peptide =  gsub(".*[.]([^.]+)[.].*", "\\1", decoy_file_2$Peptide)
  
  #combine
  combined_file = rbind(target_file_1, decoy_file_1)
  combined_file = combined_file[sample(nrow(combined_file)), ]
  combined_file = combined_file[order(-combined_file$TailorScore), ]
  
  #scan-level competition
  combined_file = combined_file[!duplicated(combined_file$ScanNr), ]
  
  #best PSM
  combined_file = combined_file[!duplicated(combined_file$Peptide), ]
  
  #get peptide-list
  peptide_list[c('decoy1', 'decoy2', 'decoy3')] = str_split_fixed(peptide_list$`decoy(s)`, ',', 3)
  
  #get original targets
  combined_file$target = combined_file$Peptide
  indxs = match(combined_file$Peptide[combined_file$Label == -1], peptide_list$decoy1)
  combined_file$target[combined_file$Label == -1] = peptide_list$target[indxs]
  
  #peptide-level competition
  combined_file = combined_file[order(-combined_file$TailorScore), ]
  combined_file = combined_file[!duplicated(combined_file$target), ]
  
  #get original targets for second decoy file
  peptide_list = peptide_list[!(peptide_list$decoy1 == peptide_list$decoy2), ]
  decoy_file_2 = decoy_file_2[decoy_file_2$Peptide %in% peptide_list$decoy2, ]
  indxs = match(decoy_file_2$Peptide, peptide_list$decoy2)
  decoy_file_2$target = peptide_list$target[indxs]
  
  #reorder decoy_file_2 to match combined_file
  combine_file_sub = combined_file[combined_file$target %in% decoy_file_2$target, ]
  indxs = match(combine_file_sub$target, decoy_file_2$target)
  decoy_file_2 = decoy_file_2[indxs, ]
  
  #create pseudolabels
  combine_file_sub$pseudolabel = 1
  decoy_file_2$pseudolabel = -1
  
  #do pseudo-competition
  SVM_train_data = decoy_file_2
  pseudo_target_wins = which(combine_file_sub$TailorScore > decoy_file_2$TailorScore)
  ties = which(combine_file_sub$TailorScore == decoy_file_2$TailorScore)
  pseudo_target_wins = c(pseudo_target_wins, ties[as.logical(sample(c(0, 1), size = length(ties), replace = TRUE))])
  SVM_train_data[pseudo_target_wins, ] = combine_file_sub[pseudo_target_wins, ]
  
  #Preprocess dataframe
  SVM_train_features = SVM_train_data %>% select(-SpecId, -Label, -ScanNr, -Peptide, -Proteins, -pseudolabel, -target)
  sds = apply(SVM_train_features, 2, sd)
  SVM_train_features = SVM_train_features[, abs(sds) > 1e-10]
  SVM_train_labels = SVM_train_data$pseudolabel
  
  #scale non-binary features
  SVM_train_features[,which(!names(SVM_train_features) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzInt", "enzN", "enzC"))] = scale(SVM_train_features[,which(!names(SVM_train_features) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzInt", "enzN", "enzC"))])
  
  #getting initial positive and negative set
  indxs = order(-SVM_train_features$TailorScore)
  SVM_train_features = SVM_train_features[indxs,]
  SVM_train_labels = SVM_train_labels[indxs]
  q_vals = TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1)
  positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= FDR_threshold)
  negative_set_indxs = (SVM_train_labels == -1)
  
  SVM_train_features_iter = SVM_train_features[positive_set_indxs | negative_set_indxs, ]
  SVM_train_labels_iter = SVM_train_labels[positive_set_indxs | negative_set_indxs]
  
  Cs_pos = Cs_neg = Cs
  
  power_finals = double()
  max_powers = double()
  
  for (iter in 1:total_iter) {
    
    #initializing C_pairs and count
    print(paste('iter:', iter))
    C_pairs = data.frame(Cs_pos = double(), Cs_neg = double())
    count = 1
    powers = double()
    #looping to find the best choice of Cs
    for (i in 1:length(Cs_neg)) {
      C_neg = Cs_neg[i]
      for (j in 1:length(Cs_pos)) {
        C_pos = Cs_pos[j]
        print('C_neg/C_pos')
        print(c(C_neg, C_pos))
        
        #train SVM
        cost_weight = c(C_pos, C_neg)
        names(cost_weight) = c(1,-1)
        
        model = svm(factor(SVM_train_labels_iter) ~.,
                    data = SVM_train_features_iter,
                    type = 'C-classification',
                    kernel = kernel,
                    class.weights = cost_weight,
                    scale = F)
        
        #determine ordering from trained SVM
        new_scores = predict(model, SVM_train_features, decision.values = TRUE)
        new_scores = as.vector(attr(new_scores, 'decision.values'))
        
        #Do TDC
        new_indxs = order(-new_scores)
        SVM_labels_new = SVM_train_labels[new_indxs]
        q_vals = TDC_flex_c(SVM_labels_new == -1, SVM_labels_new == 1)
        power = sum(q_vals <= FDR_threshold & SVM_labels_new == 1)
        
        powers = c(powers, power)
        print(paste('power', power))
        C_pairs[count, ] = cost_weight 
        
        count = count + 1
        print('yay')
      }
    }
    
    power_indx = which.max(powers)
    print('max_power')
    print(max(powers))
    max_powers = c(max_powers, max(powers))
    
    #redoing the SVM model that yielded the highest power
    #and reordering according to this model
    C_max = c(C_pairs[power_indx, 1,], C_pairs[power_indx, 2])
    names(C_max) = c(1,-1)
    
    model = svm(factor(SVM_train_labels_iter) ~.,
                data = SVM_train_features_iter,
                type = 'C-classification',
                kernel = kernel,
                class.weights = C_max,
                scale = F)
    
    new_scores = predict(model, SVM_train_features, decision.values = TRUE)
    new_scores = as.vector(attr(new_scores, 'decision.values'))
    best_indxs = order(-new_scores)
    
    SVM_train_features = SVM_train_features[best_indxs,]
    SVM_train_labels = SVM_train_labels[best_indxs]
    q_vals = TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1)
    positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= FDR_threshold)
    negative_set_indxs = (SVM_train_labels == -1)
    
    #the new positive and negative training set
    SVM_train_features_iter = SVM_train_features[positive_set_indxs | negative_set_indxs, ]
    SVM_train_labels_iter = SVM_train_labels[positive_set_indxs | negative_set_indxs]
    
    #get actual power if we were to stop here
    real_labels = combined_file$Label
    combined_file_test = combined_file %>% select(-SpecId, -Label, -ScanNr, -Peptide, -Proteins, -target)
    sds = apply(combined_file_test, 2, sd)
    combined_file_test = combined_file_test[, abs(sds) > 1e-10]
    combined_file_test[,which(!names(combined_file_test) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzInt", "enzN", "enzC"))] = scale(combined_file_test[,which(!names(combined_file_test) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzInt", "enzN", "enzC"))])
    
    new_scores = predict(model, combined_file_test, decision.values = TRUE)
    new_scores = as.vector(attr(new_scores, 'decision.values'))
    
    new_labels = real_labels[order(-new_scores)]
    
    q_val = TDC_flex_c(new_labels == -1, new_labels == 1)
    power_final = sum(q_val <= 0.01 & new_labels == 1)
    
    power_finals = c(power_finals, power_final)
  }
  
  results = list(max_powers, power_finals)
  names(results) = c('max_powers', 'power_finals')
  return(results)
}
do_iterative_svm_all = function(target_file_name_1, decoy_file_name_1, decoy_file_name_2, peptide_list_1_name, peptide_list_2_name, FDR_threshold = 0.01, Cs = c(0.1, 1, 10), total_iter = 10, kernel = 'linear') {
  
  target_file_1 = read_delim(target_file_name_1, na = character())
  decoy_file_1 = read_delim(decoy_file_name_1, na = character())
  decoy_file_2 = read_delim(decoy_file_name_2, na = character())
  peptide_list_1 = read_delim(peptide_list_1_name, na = character())
  peptide_list_2 = read_delim(peptide_list_2_name, na = character())
  
  #fix peptide format
  target_file_1$Peptide = gsub(".*[.]([^.]+)[.].*", "\\1", target_file_1$Peptide)
  decoy_file_1$Peptide = gsub(".*[.]([^.]+)[.].*", "\\1", decoy_file_1$Peptide)
  decoy_file_2$Peptide =  gsub(".*[.]([^.]+)[.].*", "\\1", decoy_file_2$Peptide)
  
  #combine
  combined_file = rbind(target_file_1, decoy_file_1)
  combined_file = combined_file[sample(nrow(combined_file)), ]
  combined_file = combined_file[order(-combined_file$TailorScore), ]
  
  #scan-level competition
  combined_file = combined_file[!duplicated(combined_file$ScanNr), ]
  
  #best PSM
  combined_file = combined_file[!duplicated(combined_file$Peptide), ]
  
  #get peptide-list
  peptide_list = peptide_list_1
  peptide_list$decoy1 = peptide_list_1$decoy
  peptide_list$decoy2 = peptide_list_2$decoy
  
  #get original targets
  combined_file$target = combined_file$Peptide
  indxs = match(combined_file$Peptide[combined_file$Label == -1], peptide_list$decoy1)
  combined_file$target[combined_file$Label == -1] = peptide_list$target[indxs]
  
  #peptide-level competition
  combined_file = combined_file[order(-combined_file$TailorScore), ]
  combined_file = combined_file[!duplicated(combined_file$target), ]
  
  #get original targets for second decoy file
  peptide_list = peptide_list[(peptide_list$decoy1 != peptide_list$decoy2) & (peptide_list$decoy1 != '') & (peptide_list$decoy2 != ''), ]
  decoy_file_2 = decoy_file_2[decoy_file_2$Peptide %in% peptide_list$decoy2, ]
  indxs = match(decoy_file_2$Peptide, peptide_list$decoy2)
  decoy_file_2$target = peptide_list$target[indxs]
  
  #make target column
  combine_file_sub = combined_file[combined_file$target %in% peptide_list$target, ]
  combine_file_sub$target = combine_file_sub$Peptide

  #create pseudolabels
  combine_file_sub$pseudolabel = 1
  decoy_file_2$pseudolabel = -1
  
  SVM_train_data = rbind(combine_file_sub, decoy_file_2)
  SVM_train_data = SVM_train_data[sample(nrow(SVM_train_data)), ]
  SVM_train_data = SVM_train_data[order(-SVM_train_data$TailorScore), ]
  SVM_train_data = SVM_train_data[!duplicated(SVM_train_data$target), ]
  
  #Preprocess dataframe
  SVM_train_features = SVM_train_data %>% select(-SpecId, -Label, -ScanNr, -Peptide, -Proteins, -pseudolabel, -target)
  sds = apply(SVM_train_features, 2, sd)
  SVM_train_features = SVM_train_features[, abs(sds) > 1e-10]
  SVM_train_labels = SVM_train_data$pseudolabel
  
  #scale non-binary features
  SVM_train_features[,which(!names(SVM_train_features) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzInt", "enzN", "enzC"))] = scale(SVM_train_features[,which(!names(SVM_train_features) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzInt", "enzN", "enzC"))])

  #getting initial positive and negative set
  indxs = order(-SVM_train_features$TailorScore)
  SVM_train_features = SVM_train_features[indxs,]
  SVM_train_labels = SVM_train_labels[indxs]
  q_vals = TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1)
  positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= FDR_threshold)
  negative_set_indxs = (SVM_train_labels == -1)
  
  SVM_train_features_iter = SVM_train_features[positive_set_indxs | negative_set_indxs, ]
  SVM_train_labels_iter = SVM_train_labels[positive_set_indxs | negative_set_indxs]
  
  Cs_pos = Cs_neg = Cs
  
  power_finals = double()
  max_powers = double()
  
  for (iter in 1:total_iter) {
    
    #initializing C_pairs and count
    print(paste('iter:', iter))
    C_pairs = data.frame(Cs_pos = double(), Cs_neg = double())
    count = 1
    powers = double()
    #looping to find the best choice of Cs
    for (i in 1:length(Cs_neg)) {
      C_neg = Cs_neg[i]
      for (j in 1:length(Cs_pos)) {
        C_pos = Cs_pos[j]
        print('C_neg/C_pos')
        print(c(C_neg, C_pos))
        
        #train SVM
        cost_weight = c(C_pos, C_neg)
        names(cost_weight) = c(1,-1)
        
        model = svm(factor(SVM_train_labels_iter) ~.,
                    data = SVM_train_features_iter,
                    type = 'C-classification',
                    kernel = kernel,
                    class.weights = cost_weight,
                    scale = T)
        
        #determine ordering from trained SVM
        new_scores = predict(model, SVM_train_features, decision.values = TRUE)
        new_scores = as.vector(attr(new_scores, 'decision.values'))
        
        #Do TDC
        new_indxs = order(-new_scores)
        SVM_labels_new = SVM_train_labels[new_indxs]
        q_vals = TDC_flex_c(SVM_labels_new == -1, SVM_labels_new == 1)
        power = sum(q_vals <= FDR_threshold & SVM_labels_new == 1)
        
        powers = c(powers, power)
        print(paste('power', power))
        C_pairs[count, ] = cost_weight 
  
        count = count + 1
        print('yay')
      }
    }
    
    power_indx = which.max(powers)
    print('max_power')
    print(max(powers))
    max_powers = c(max_powers, max(powers))
    
    #redoing the SVM model that yielded the highest power
    #and reordering according to this model
    C_max = c(C_pairs[power_indx, 1,], C_pairs[power_indx, 2])
    names(C_max) = c(1,-1)
    
    model = svm(factor(SVM_train_labels_iter) ~.,
                data = SVM_train_features_iter,
                type = 'C-classification',
                kernel = kernel,
                class.weights = C_max,
                scale = T)
    
    new_scores = predict(model, SVM_train_features, decision.values = TRUE)
    new_scores = as.vector(attr(new_scores, 'decision.values'))
    best_indxs = order(-new_scores)
    
    SVM_train_features = SVM_train_features[best_indxs,]
    SVM_train_labels = SVM_train_labels[best_indxs]
    q_vals = TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1)
    positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= FDR_threshold)
    negative_set_indxs = (SVM_train_labels == -1)
    
    #the new positive and negative training set
    SVM_train_features_iter = SVM_train_features[positive_set_indxs | negative_set_indxs, ]
    SVM_train_labels_iter = SVM_train_labels[positive_set_indxs | negative_set_indxs]
    
    #get actual power if we were to stop here
    real_labels = combined_file$Label
    combined_file_test = combined_file %>% select(-SpecId, -Label, -ScanNr, -Peptide, -Proteins, -target)
    sds = apply(combined_file_test, 2, sd)
    combined_file_test = combined_file_test[, abs(sds) > 1e-10]
    combined_file_test[,which(!names(combined_file_test) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzInt", "enzN", "enzC"))] = scale(combined_file_test[,which(!names(combined_file_test) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzInt", "enzN", "enzC"))])
    
    new_scores = predict(model, combined_file_test, decision.values = TRUE)
    new_scores = as.vector(attr(new_scores, 'decision.values'))
    
    new_labels = real_labels[order(-new_scores)]
    
    q_val = TDC_flex_c(new_labels == -1, new_labels == 1)
    power_final = sum(q_val <= 0.01 & new_labels == 1)
    
    power_finals = c(power_finals, power_final)
    print('POWER_FINAL')
  }
  
  results = list(max_powers, power_finals)
  names(results) = c('max_powers', 'power_finals')
  return(results)
}
do_iterative_svm_dd_all = function(target_file_name_1, decoy_file_name_1, decoy_file_name_2, decoy_file_name_3, peptide_list_1_name, peptide_list_2_name, peptide_list_3_name, FDR_threshold = 0.01, Cs = c(0.1, 1, 10), total_iter = 10, kernel = 'linear') {
  
  #reading in the files
  target_file_1 = read_delim(target_file_name_1, na = character())
  decoy_file_1 = read_delim(decoy_file_name_1, na = character())
  decoy_file_2 = read_delim(decoy_file_name_2, na = character())
  decoy_file_3 = read_delim(decoy_file_name_3, na = character())
  peptide_list_1 = read_delim(peptide_list_name_1, na = character())
  peptide_list_2 = read_delim(peptide_list_name_2, na = character())
  peptide_list_3 = read_delim(peptide_list_name_3, na = character())
  
  #fix peptide format
  target_file_1$Peptide = gsub(".*[.]([^.]+)[.].*", "\\1", target_file_1$Peptide)
  decoy_file_1$Peptide = gsub(".*[.]([^.]+)[.].*", "\\1", decoy_file_1$Peptide)
  decoy_file_2$Peptide =  gsub(".*[.]([^.]+)[.].*", "\\1", decoy_file_2$Peptide)
  decoy_file_3$Peptide =  gsub(".*[.]([^.]+)[.].*", "\\1", decoy_file_3$Peptide)
  
  #combine
  combined_file = rbind(target_file_1, decoy_file_1)
  combined_file = combined_file[sample(nrow(combined_file)), ]
  combined_file = combined_file[order(-combined_file$TailorScore), ]
  
  #scan-level competition
  combined_file = combined_file[!duplicated(combined_file$ScanNr), ]
  
  #best PSM
  combined_file = combined_file[!duplicated(combined_file$Peptide), ]
  
  #get peptide-list
  peptide_list[c('decoy1', 'decoy2', 'decoy3')] = str_split_fixed(peptide_list$`decoy(s)`, ',', 3)
  
  #get original targets
  combined_file$target = combined_file$Peptide
  indxs = match(combined_file$Peptide[combined_file$Label == -1], peptide_list$decoy1)
  combined_file$target[combined_file$Label == -1] = peptide_list$target[indxs]
  
  #peptide-level competition
  combined_file = combined_file[order(-combined_file$TailorScore), ]
  combined_file = combined_file[!duplicated(combined_file$target), ]
  
  #use only peptides with distinct decoys
  peptide_list = peptide_list[str_count(peptide_list$`decoy(s)`, ',') == 2, ]
  
  #do spectrum-level competition for decoy_file_2/3
  decoy_file_2$decoy = 'decoy2'
  decoy_file_3$decoy = 'decoy3'
  decoy_file_2 = decoy_file_2[decoy_file_2$Peptide %in% peptide_list$decoy2,]
  decoy_file_3 = decoy_file_3[decoy_file_3$Peptide %in% peptide_list$decoy3,]
  decoy_file_combined = rbind(decoy_file_2, decoy_file_3)
  decoy_file_combined = decoy_file_combined[sample(nrow(decoy_file_combined)), ]
  decoy_file_combined = decoy_file_combined[order(-decoy_file_combined$TailorScore), ]
  decoy_file_combined = decoy_file_combined[!duplicated(decoy_file_combined$ScanNr), ]
  decoy_file_combined = decoy_file_combined[!duplicated(decoy_file_combined$Peptide), ]
  
  #get original targets for decoy_file_combined
  decoy_file_combined$target = decoy_file_combined$Peptide
  indxs = match(decoy_file_combined$Peptide[decoy_file_combined$decoy == 'decoy2'], peptide_list$decoy2)
  decoy_file_combined$target[decoy_file_combined$decoy == 'decoy2'] = peptide_list$target[indxs]
  indxs = match(decoy_file_combined$Peptide[decoy_file_combined$decoy == 'decoy3'], peptide_list$decoy3)
  decoy_file_combined$target[decoy_file_combined$decoy == 'decoy3'] = peptide_list$target[indxs]
  
  #do peptide-level competition
  decoy_file_combined = decoy_file_combined[order(-decoy_file_combined$TailorScore), ]
  decoy_file_combined = decoy_file_combined[!duplicated(decoy_file_combined$target), ]
  
  #reorder decoy_file_combined to match combined_file
  combine_file_sub = combined_file[combined_file$target %in% decoy_file_combined$target, ]
  indxs = match(combine_file_sub$target, decoy_file_combined$target)
  decoy_file_combined = decoy_file_combined[indxs, ]
  decoy_file_combined = decoy_file_combined %>% select(-decoy)
  
  #create pseudolabels
  combine_file_sub$pseudolabel = 1
  decoy_file_combined$pseudolabel = -1
  
  #do pseudo-competition
  SVM_train_data = decoy_file_combined
  pseudo_target_wins = which(combine_file_sub$TailorScore > decoy_file_combined$TailorScore)
  ties = which(combine_file_sub$TailorScore == decoy_file_combined$TailorScore)
  pseudo_target_wins = c(pseudo_target_wins, ties[as.logical(sample(c(0, 1), size = length(ties), replace = TRUE))])
  SVM_train_data[pseudo_target_wins, ] = combine_file_sub[pseudo_target_wins, ]
  
  #Preprocess dataframe
  SVM_train_features = SVM_train_data %>% select(-SpecId, -Label, -ScanNr, -Peptide, -Proteins, -pseudolabel, -target)
  sds = apply(SVM_train_features, 2, sd)
  SVM_train_features = SVM_train_features[, abs(sds) > 1e-10]
  SVM_train_labels = SVM_train_data$pseudolabel
  
  #scale non-binary features
  SVM_train_features[,which(!names(SVM_train_features) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzInt", "enzN", "enzC"))] = scale(SVM_train_features[,which(!names(SVM_train_features) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzInt", "enzN", "enzC"))])
  
  #getting initial positive and negative set
  indxs = order(-SVM_train_features$TailorScore)
  SVM_train_features = SVM_train_features[indxs,]
  SVM_train_labels = SVM_train_labels[indxs]
  q_vals = TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1)
  positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= FDR_threshold)
  negative_set_indxs = (SVM_train_labels == -1)
  
  SVM_train_features_iter = SVM_train_features[positive_set_indxs | negative_set_indxs, ]
  SVM_train_labels_iter = SVM_train_labels[positive_set_indxs | negative_set_indxs]
  
  C_pos = C_neg = Cs
  
  power_finals = double()
  max_powers = double()
  
  for (iter in 1:total_iter) {
    
    #initializing C_pairs and count
    print(paste('iter:', iter))
    C_pairs = data.frame(Cs_pos = double(), Cs_neg = double())
    count = 1
    powers = double()
    #looping to find the best choice of Cs
    for (i in 1:length(Cs_neg)) {
      C_neg = Cs_neg[i]
      for (j in 1:length(Cs_pos)) {
        C_pos = Cs_pos[j]
        print('C_neg/C_pos')
        print(c(C_neg, C_pos))
        
        #train SVM
        cost_weight = c(C_pos, C_neg)
        names(cost_weight) = c(1,-1)
        
        model = svm(factor(SVM_train_labels_iter) ~.,
                    data = SVM_train_features_iter,
                    type = 'C-classification',
                    kernel = kernel,
                    class.weights = cost_weight,
                    scale = F)
        
        #determine ordering from trained SVM
        new_scores = predict(model, SVM_train_features, decision.values = TRUE)
        new_scores = as.vector(attr(new_scores, 'decision.values'))
        
        #Do TDC
        new_indxs = order(-new_scores)
        SVM_labels_new = SVM_train_labels[new_indxs]
        q_vals = TDC_flex_c(SVM_labels_new == -1, SVM_labels_new == 1)
        power = sum(q_vals <= FDR_threshold & SVM_labels_new == 1)
        
        powers = c(powers, power)
        print(paste('power', power))
        C_pairs[count, ] = cost_weight 
        
        count = count + 1
        print('yay')
      }
    }
    
    power_indx = which.max(powers)
    print('max_power')
    print(max(powers))
    max_powers = c(max_powers, max(powers))
    
    #redoing the SVM model that yielded the highest power
    #and reordering according to this model
    C_max = c(C_pairs[power_indx, 1,], C_pairs[power_indx, 2])
    names(C_max) = c(1,-1)
    
    model = svm(factor(SVM_train_labels_iter) ~.,
                data = SVM_train_features_iter,
                type = 'C-classification',
                kernel = kernel,
                class.weights = C_max,
                scale = F)
    
    new_scores = predict(model, SVM_train_features, decision.values = TRUE)
    new_scores = as.vector(attr(new_scores, 'decision.values'))
    best_indxs = order(-new_scores)
    
    SVM_train_features = SVM_train_features[best_indxs,]
    SVM_train_labels = SVM_train_labels[best_indxs]
    q_vals = TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1)
    positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= FDR_threshold)
    negative_set_indxs = (SVM_train_labels == -1)
    
    #the new positive and negative training set
    SVM_train_features_iter = SVM_train_features[positive_set_indxs | negative_set_indxs, ]
    SVM_train_labels_iter = SVM_train_labels[positive_set_indxs | negative_set_indxs]
    
    #get actual power if we were to stop here
    real_labels = combined_file$Label
    combined_file_test = combined_file %>% select(-SpecId, -Label, -ScanNr, -Peptide, -Proteins, -target)
    sds = apply(combined_file_test, 2, sd)
    combined_file_test = combined_file_test[, abs(sds) > 1e-10]
    combined_file_test[,which(!names(combined_file_test) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzInt", "enzN", "enzC"))] = scale(combined_file_test[,which(!names(combined_file_test) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzInt", "enzN", "enzC"))])
    
    new_scores = predict(model, combined_file_test, decision.values = TRUE)
    new_scores = as.vector(attr(new_scores, 'decision.values'))
    
    new_labels = real_labels[order(-new_scores)]
    
    q_val = TDC_flex_c(new_labels == -1, new_labels == 1)
    power_final = sum(q_val <= 0.01 & new_labels == 1)
    
    power_finals = c(power_finals, power_final)
  }
  
  results = list(max_powers, power_finals)
  names(results) = c('max_powers', 'power_finals')
  return(results)
}

#Procedure in a nutshell
# 1. Start with estimating the top 1% discoveries 
# 2. Re-estimate the new direction
# 3. Repeats steps 1 and 2, except starting with the new direction

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
peptide_level = function(narrow_file, open_file, peptide_list) { #just for static modifications
  narrow_df = read_delim(narrow_file, na = character(), show_col_types = FALSE)
  open_df = read_delim(open_file, na = character(), show_col_types = FALSE)
  narrow_df$n_o = 1
  open_df$n_o = 0
  
  df_all = rbind(narrow_df, open_df)
  df_all = df_all[sample(nrow(df_all)), ]
  df_all = df_all[order(-df_all$TailorScore), ]
  df_all$rank = unlist(lapply(df_all$SpecId, function(x) substr(x, start = nchar(x), nchar(x))))
  df_all = df_all[df_all$rank <= 2, ]
  df_all$rank[df_all$rank == 1] = 0 #transforming into a binary
  df_all$rank[df_all$rank == 2] = 1
  df_all$rank = as.numeric(df_all$rank)
  
  peptide_list_df = read_delim(peptide_list, na = character())
  original_target = df_all$Peptide
  original_target = unlist(lapply(original_target, function(x) substr(x, start = 3, nchar(x) - 2)))
  original_target[df_all$Label == -1] = peptide_list_df$target[match(original_target[df_all$Label == -1], peptide_list_df$decoy)]
  
  df_all = df_all[!duplicated(original_target), ]
  
  df_all = df_all %>% select(-enzInt)
  return(df_all)
}
do_iterative_svm = function(df_all, train_prob = 0.5, Cs = c(0.1, 1, 10), total_iter = 5, kernel = 'linear', alpha = 0.01) {
  #create train dataframe
  train_decoys_indx = sample(c(T, F), size = sum(df_all$Label == -1), replace = T, prob = c(train_prob, 1 - train_prob))
  train_decoys = df_all[which(df_all$Label == -1)[train_decoys_indx], ]
  train_targets = df_all[-which(df_all$Label == -1)[train_decoys_indx], ]
  train_targets$Label = 1
  train_df = rbind(train_decoys, train_targets)
  
  #real dataframe
  real_df = df_all[-which(df_all$Label == -1)[train_decoys_indx], ]
  
  #Preprocess dataframe
  SVM_train_labels = train_df$Label
  SVM_train_features = train_df %>% select(-SpecId, -Label, -ScanNr, -Peptide, -Proteins)
  sds = apply(SVM_train_features, 2, sd)
  SVM_train_features = SVM_train_features[, abs(sds) > 1e-10]
  
  #scale non-binary features
  SVM_train_features[,which(!names(SVM_train_features) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzN", "enzC", "rank"))] = scale(SVM_train_features[,which(!names(SVM_train_features) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzN", "enzC", "rank"))])
  
  #getting initial positive and negative set
  indxs = order(-SVM_train_features$TailorScore)
  SVM_train_features = SVM_train_features[indxs,]
  SVM_train_labels = SVM_train_labels[indxs]
  q_vals = TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1, c = 3/4, lambda = 3/4)
  positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= alpha)
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
                    degree = degree,
                    class.weights = cost_weight,
                    scale = F)
        
        #determine ordering from trained SVM
        new_scores = predict(model, SVM_train_features, decision.values = TRUE)
        new_scores = as.vector(attr(new_scores, 'decision.values'))
        
        #Do TDC
        new_indxs = order(-new_scores)
        SVM_labels_new = SVM_train_labels[new_indxs]
        q_vals = TDC_flex_c(SVM_labels_new == -1, SVM_labels_new == 1, c = 3/4, lambda = 3/4)
        power = sum(q_vals <= alpha & SVM_labels_new == 1)
        
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
    q_vals = TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1, c = 3/4, lambda = 3/4)
    positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= alpha)
    negative_set_indxs = (SVM_train_labels == -1)
    
    #the new positive and negative training set
    SVM_train_features_iter = SVM_train_features[positive_set_indxs | negative_set_indxs, ]
    SVM_train_labels_iter = SVM_train_labels[positive_set_indxs | negative_set_indxs]
    
    #get actual power if we were to stop here
    real_labels = real_df$Label
    real_df_test = real_df %>% select(-SpecId, -Label, -ScanNr, -Peptide, -Proteins)
    sds = apply(real_df_test, 2, sd)
    real_df_test = real_df_test[, abs(sds) > 1e-10]
    real_df_test[,which(!names(real_df_test) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzN", "enzC", "rank"))] = scale(real_df_test[,which(!names(real_df_test) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzN", "enzC", "rank"))])
    
    new_scores = predict(model, real_df_test, decision.values = TRUE)
    new_scores = as.vector(attr(new_scores, 'decision.values'))
    
    new_labels = real_labels[order(-new_scores)]
    
    q_val = TDC_flex_c(new_labels == -1, new_labels == 1, c = 2/3, lambda = 2/3)
    power_final = sum(q_val <= 0.01 & new_labels == 1)
    
    power_finals = c(power_finals, power_final)
  }
  results = list(max_powers, power_finals)
  names(results) = c('max_powers', 'power_finals')
  return(results)
}
train_cv = function(labels, df, folds = 3, Cs = c(0.1, 1, 10), kernel = 'linear', degree = NULL, alpha = 0.01) {
  Cs_neg = Cs_pos = Cs
  
  indx_sample = sample(0:folds, size = nrow(df), replace = T)
  
  df_results = data.frame(
    C_pos = double(),
    C_neg = double(),
    average_power = double()
  )
  
  count_df = 1
  for (j in 1:length(Cs_pos)) {
    for (i in j:length(Cs_neg)) { #negative should be at least as large as positive weight? Penalizing negative is more important?
      C_neg = Cs_neg[i]
      C_pos = Cs_pos[j]
      print('C_neg/C_pos')
      print(c(C_neg, C_pos))
      
      #train SVM
      cost_weight = c(C_pos, C_neg)
      names(cost_weight) = c(1,-1)
      
      power_tests = double(length = folds)
      
      count = 1
      for (fold in 0:folds) {
        train_folds = setdiff(0:folds, fold)
        model = svm(factor(labels[indx_sample %in% train_folds]) ~.,
                    data = df[indx_sample %in% train_folds, ],
                    type = 'C-classification',
                    kernel = kernel,
                    degree = degree,
                    class.weights = cost_weight,
                    scale = F)
        
        #determine ordering from trained SVM
        new_scores = predict(model, df[indx_sample == fold, ], decision.values = TRUE)
        new_scores = as.vector(attr(new_scores, 'decision.values'))
        
        #Do TDC
        new_indxs = order(-new_scores)
        SVM_labels_new = labels[indx_sample == fold][new_indxs]
        q_vals = TDC_flex_c(SVM_labels_new == -1, SVM_labels_new == 1, c = 3/4, lambda = 3/4)
        power = sum(q_vals <= alpha & SVM_labels_new == 1)
        
        power_tests[count] = power
        count = count + 1
      }
      
      df_results[count_df, ] = c(C_pos, C_neg, mean(power_tests))
    
      count_df = count_df + 1
    }
  }
  return(df_results)
}
do_iterative_svm_cv = function(df_all, folds = 3, train_prob = 0.5, Cs = c(0.1, 1, 10), total_iter = 10, kernel = 'linear', alpha = 0.01, train_alpha = 0.01, degree = NULL) {
  #create train dataframe
  train_decoys_indx = sample(c(T, F), size = sum(df_all$Label == -1), replace = T, prob = c(train_prob, 1 - train_prob))
  train_decoys = df_all[which(df_all$Label == -1)[train_decoys_indx], ]
  train_targets = df_all[-which(df_all$Label == -1)[train_decoys_indx], ]
  train_targets$Label = 1
  train_df = rbind(train_decoys, train_targets)
  
  #real dataframe
  real_df = df_all[-which(df_all$Label == -1)[train_decoys_indx], ]
  
  #Preprocess dataframe
  SVM_train_labels = train_df$Label
  SVM_train_features = train_df %>% select(-SpecId, -Label, -ScanNr, -Peptide, -Proteins)
  sds = apply(SVM_train_features, 2, sd)
  SVM_train_features = SVM_train_features[, abs(sds) > 1e-10]
  
  #scale non-binary features
  SVM_train_features[,which(!names(SVM_train_features) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzN", "enzC", "rank", "n_o"))] = scale(SVM_train_features[,which(!names(SVM_train_features) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzN", "enzC", "rank", "n_o"))])
  
  #getting initial positive and negative set
  indxs = order(-SVM_train_features$TailorScore)
  SVM_train_features = SVM_train_features[indxs,]
  SVM_train_labels = SVM_train_labels[indxs]
  q_vals = TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1, c = 3/4, lambda = 3/4)
  positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= alpha)
  negative_set_indxs = (SVM_train_labels == -1)
  
  SVM_train_features_iter = SVM_train_features[positive_set_indxs | negative_set_indxs, ]
  SVM_train_labels_iter = SVM_train_labels[positive_set_indxs | negative_set_indxs]
  
  Cs_pos = Cs_neg = Cs
  
  power_finals = double()
  max_powers = double()
  
  for (iter in 1:total_iter) {
   
    print(paste('iter:', iter))
    
    results_df = train_cv(SVM_train_labels_iter, SVM_train_features_iter, folds = folds, Cs = Cs, kernel = kernel, degree = degree, alpha = train_alpha)
    
    power_indx = which.max(results_df$average_power)
    print('max_power')
    print(max(results_df$average_power))
    max_powers = c(max_powers, max(results_df$average_power))
    
    
    #redoing the SVM model that yielded the highest power
    #and reordering according to this model
    C_max = c(results_df$C_pos[power_indx], results_df$C_neg[power_indx])
    names(C_max) = c(1,-1)
    
    model = svm(factor(SVM_train_labels_iter) ~.,
                data = SVM_train_features_iter,
                type = 'C-classification',
                kernel = kernel,
                class.weights = C_max,
                degree = degree,
                scale = F)
    
    new_scores = predict(model, SVM_train_features, decision.values = TRUE)
    new_scores = as.vector(attr(new_scores, 'decision.values'))
    best_indxs = order(-new_scores)
    
    SVM_train_features = SVM_train_features[best_indxs,]
    SVM_train_labels = SVM_train_labels[best_indxs]
    q_vals = TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1, c = 3/4, lambda = 3/4)
    positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= alpha)
    negative_set_indxs = (SVM_train_labels == -1)
    
    #the new positive and negative training set
    SVM_train_features_iter = SVM_train_features[positive_set_indxs | negative_set_indxs, ]
    SVM_train_labels_iter = SVM_train_labels[positive_set_indxs | negative_set_indxs]
    
    #get actual power if we were to stop here
    real_labels = real_df$Label
    real_df_test = real_df %>% select(-SpecId, -Label, -ScanNr, -Peptide, -Proteins)
    sds = apply(real_df_test, 2, sd)
    real_df_test = real_df_test[, abs(sds) > 1e-10]
    real_df_test[,which(!names(real_df_test) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzN", "enzC", "rank", "n_o"))] = scale(real_df_test[,which(!names(real_df_test) %in% c("Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzN", "enzC", "rank", "n_o"))])
    
    new_scores = predict(model, real_df_test, decision.values = TRUE)
    new_scores = as.vector(attr(new_scores, 'decision.values'))
    
    new_labels = real_labels[order(-new_scores)]
    
    q_val = TDC_flex_c(new_labels == -1, new_labels == 1, c = 2/3, lambda = 2/3)
    power_final = sum(q_val <= 0.01 & new_labels == 1)
    
    power_finals = c(power_finals, power_final)
    print(paste("power_final:", power_final))
  }
  results = list(max_powers, power_finals)
  names(results) = c('max_powers', 'power_finals')
  return(results)
}


narrow_file = 'narrow_filtered.make-pin.pin'
open_file = 'open_filtered.make-pin.pin'
peptide_list = 'tide-index.peptides.txt'

set.seed(1)
#PXD025130
#CONGA gives 7174 discoveries
#Percolator gives 6949 discoveries
#Super Percolator gives 7219 discoveries

df_all = peptide_level(narrow_file, open_file, peptide_list)
results = do_iterative_svm_cv(df_all, train_alpha = 0.01)

set.seed(1)
#PXD025130
#CONGA gives 7174 discoveries
#Percolator gives 6949 discoveries
#Super Percolator gives 7245 discoveries with train_alpha = 0.008

df_all = peptide_level(narrow_file, open_file, peptide_list)
results = do_iterative_svm_cv(df_all, train_alpha = 0.008) #since there is a larger variability in the train data, we lower the FDR threshold

set.seed(1)
#PXD025130
#CONGA gives 7174 discoveries
#Percolator gives 6949 discoveries
#Super Percolator gives 7236 discoveries with train_alpha = 0.005

df_all = peptide_level(narrow_file, open_file, peptide_list)
results = do_iterative_svm_cv(df_all, train_alpha = 0.01, kernel = 'polynomial', degree = 2) #since there is a larger variability in the train data, we lower the FDR threshold

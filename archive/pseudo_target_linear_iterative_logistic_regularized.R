#setting wd
setwd("~/Documents/PhD/ppx_files/High resolution/Static_modification/PXD025130/crux-output/super_percolator")

#libraries
library(tidyverse)
library(caret)
library(e1071)
library(CVXR)

#parameters
FDR_threshold = 0.01

#functions
TDC_flex_c = function(decoy_wins, target_wins, BC1 = 1, c = 1/2, lambda = 1/2){
  #give in order of best score first to worst score last
  nTD <- cumsum(target_wins)
  nDD <- cumsum(decoy_wins)
  fdps <- (pmin(1, ((BC1 + nDD)/ nTD) * (c / (1-lambda)))) #nDD, nTD have same length
  qvals <- rev(cummin(rev(fdps)))
  return(qvals)
}

#read in percolator file from search-0
target_file = read_delim('./search-0/enzint_zero.tide-search.target.pin', delim = '\t')
decoy_file = read_delim('./search-0/enzint_zero.tide-search.decoy.pin', delim = '\t')

#read in decoy percolator file from search-1
decoy_file_2 = read_delim('./search-1/enzint_zero.tide-search.decoy.pin', delim = '\t')

target_file$Peptide = gsub(".*[.]([^.]+)[.].*", "\\1", target_file$Peptide)
decoy_file$Peptide = gsub(".*[.]([^.]+)[.].*", "\\1", decoy_file$Peptide)
decoy_file_2$Peptide =  gsub(".*[.]([^.]+)[.].*", "\\1", decoy_file_2$Peptide)

#combining
combined_file = rbind(target_file, decoy_file)
combined_file = combined_file[sample(nrow(combined_file)), ]
combined_file = combined_file[order(-combined_file$TailorScore), ]
combined_file = combined_file[!duplicated(combined_file$ScanNr), ]
combined_file = combined_file[!duplicated(combined_file$Peptide), ]

#reading in peptide-list
peptide_list_0 = read.delim('./index-0/tide-index.peptides.txt')
peptide_list_1 = read.delim('./index-1/tide-index.peptides.txt')

peptide_list = data.frame(peptide_list_0$target, peptide_list_0$decoy, peptide_list_1$decoy)
names(peptide_list) = c('target', 'decoy1', 'decoy2')

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
#not that many pairs, might need to go to more than top 1 rank
combine_file_sub = combined_file[combined_file$target %in% decoy_file_2$target, ]
indxs = match(combine_file_sub$target, decoy_file_2$target)
decoy_file_2 = decoy_file_2[indxs, ]

#create pseudolabels
combine_file_sub$pseudolabel = 1
decoy_file_2$pseudolabel = -1

#do pseudo-competition
LR_train_data = decoy_file_2
pseudo_target_wins = which(combine_file_sub$TailorScore > decoy_file_2$TailorScore)
ties = which(combine_file_sub$TailorScore == decoy_file_2$TailorScore)
pseudo_target_wins = c(pseudo_target_wins, ties[as.logical(sample(c(0, 1), size = length(ties), replace = TRUE))])
LR_train_data[pseudo_target_wins, ] = combine_file_sub[pseudo_target_wins, ]

#Preprocess dataframe
LR_train_features = LR_train_data %>% select(pseudolabel, ExpMass, CalcMass, deltCn, XCorr, TailorScore,
                                               PepLen, Charge2, Charge3, lnNumSP, dM, absdM,
)
LR_labels = LR_train_data$pseudolabel

attach(LR_train_features)
LR_train_features = LR_train_features %>% select(-Charge3) %>% scale() %>% as.data.frame()

total_iter = 5

indxs = order(-LR_train_features[, 5])
LR_train_features = LR_train_features[indxs,]
LR_labels = LR_labels[indxs]
q_vals = TDC_flex_c(LR_labels == -1, LR_labels == 1)
positive_set_indxs = (LR_labels == 1) & (q_vals <= FDR_threshold)
negative_set_indxs = (LR_labels == -1)

LR_train_features_iter = LR_train_features[positive_set_indxs | negative_set_indxs, ]

x = model.matrix(pseudolabel~., LR_train_features_iter)[,-1]
y = ifelse(LR_train_features_iter$pseudolabel > 0, 1, 0)

power_finals = double()
power_maxs = double()

for (iter in 1:10) {
  print(paste('iter:', iter))
  
  lambdas = glmnet(x, y, alpha = 0, family = "binomial")$lambda
  powers = rep(0, 100)
  count = 1
  for (lambda in lambdas){
    model = glmnet(x, y, alpha = 0, family = "binomial", lambda = lambda)
    
    x.test = model.matrix(pseudolabel~., LR_train_features)[,-1]
    label.test = LR_train_features$pseudolabel
    new_scores = model %>% predict(x.test)
      
    label.test = label.test[order(-new_scores)]
    
    q_vals = TDC_flex_c(label.test < 0, label.test > 0)
    powers[count] = sum((label.test > 0) & (q_vals <= FDR_threshold))
    count = count + 1
  }
  power_maxs = c(power_maxs, max(powers))
  power_indx = which.max(powers)
  model = glmnet(x, y, alpha = 0, family = "binomial", lambda = lambdas[power_indx])
  
  x.test = model.matrix(pseudolabel~., LR_train_features)[,-1]
  new_scores = model %>% predict(x.test)
  
  LR_train_features = LR_train_features[order(-new_scores), ]
  LR_labels = LR_train_features$pseudolabel
  
  q_vals = TDC_flex_c(LR_labels < 0, LR_labels > 0)
  positive_set_indxs = (LR_labels > 0) & (q_vals <= FDR_threshold)
  negative_set_indxs = (LR_labels < 0)
  
  LR_train_features_iter = LR_train_features[positive_set_indxs | negative_set_indxs, ]
  
  x = model.matrix(pseudolabel~., LR_train_features_iter)[,-1]
  y = ifelse(LR_train_features_iter$pseudolabel > 0, 1, 0)
  
  ######################################################################################
  
  combined_file_features = combined_file %>% select(Label, ExpMass, CalcMass, deltCn, XCorr, TailorScore,
                                                    PepLen, Charge2, Charge3, lnNumSP, dM, absdM,
  )
  attach(combined_file_features)
  combined_file_features = combined_file_features %>% select(-Charge3) %>% scale() %>% as.data.frame()
  
  x.final = model.matrix(Label~., combined_file_features)[, -1]
  Label.final = combined_file$Label
  new_scores = model %>% predict(x.final)
  Label.final = Label.final[order(-new_scores)]
  
  q_val = TDC_flex_c(Label.final == -1, Label.final == 1)
  power_final = sum(q_val <= 0.01 & Label.final == 1)
  power_finals = c(power_finals, power_final)
  
  print(power_final)
  Sys.sleep(1)
}

#compare this to narrow TDC
new_labels = combined_file$Label[order(-combined_file$TailorScore)]
q_val = TDC_flex_c(new_labels == -1, new_labels == 1)
power_final = sum(q_val <= 0.01 & new_labels == 1) #4712 discoveries

#percolator reports 4821


#setting wd
setwd("~/Documents/PhD/ppx_files/High resolution/Static_modification/PXD025130/crux-output/super_percolator")

#libraries
library(tidyverse)
library(caret)
library(e1071)
library(CVXR)

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
combined_file = combined_file[order(-combined_file$TailorScore), ]
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
SVM_train_data = decoy_file_2
pseudo_target_wins = which(combine_file_sub$TailorScore > decoy_file_2$TailorScore)
ties = which(combine_file_sub$TailorScore == decoy_file_2$TailorScore)
pseudo_target_wins = c(pseudo_target_wins, ties[as.logical(sample(c(0, 1), size = length(ties), replace = TRUE))])
SVM_train_data[pseudo_target_wins, ] = combine_file_sub[pseudo_target_wins, ]

#Preprocess dataframe
SVM_train_features = SVM_train_data %>% select(ExpMass, CalcMass, deltCn, XCorr, TailorScore,
                                               PepLen, Charge2, Charge3, lnNumSP, dM, absdM,
                                               )
SVM_labels = SVM_train_data$pseudolabel

attach(SVM_train_features)
SVM_train_features = SVM_train_features %>% select(-Charge2, -Charge3) %>% scale() %>% as.data.frame()
SVM_train_features = cbind(SVM_train_features, Charge2) 
SVM_train_features = as.matrix(SVM_train_features)

# If you want to take a subset for testing, use this
# inds = sample(nrow(SVM_train_features), size = 1000, replace = F)
# SVM_train_features = SVM_train_features[inds, ]
# SVM_labels = SVM_labels[inds]


####################################################

Cs_pos = c(0.01, 0.1, 1, 10, 100, 1000)
Cs_neg = c(0.01, 0.1, 1, 10, 100, 1000)

coef_df = data.frame(Cs_pos = double(), Cs_neg = double(), 
                     betas = replicate(ncol(SVM_train_features) + 1, double())
                     )

count = 1
for (i in 1:length(Cs_neg)) {
  C_neg = Cs_neg[i]
  for (j in 1:i) {
    C_pos = Cs_pos[j]
    print(c(C_neg, C_pos))
    n = ncol(SVM_train_features)
    m = nrow(SVM_train_features)
    C = rep(0, m)
    C[SVM_labels == 1] = C_pos
    C[SVM_labels == -1] = C_neg
    C = matrix(C, nrow = 1)
    
    beta = Variable(n)
    beta_0 = Variable(1)

    obj = C%*%pos(1 - SVM_labels * (SVM_train_features%*%beta + beta_0 )) + 0.5 * sum(beta^2)
    prob = Problem(Minimize(obj))
    result = solve(prob, verbose = T, solver = 'ECOS') #ECOS appears to be the most efficient and reliable solver
    
    coef_df[count, ] = c(C_pos, C_neg, result$getValue(beta), result$getValue(beta_0))
    count = count + 1
    print('yay')
  }
}

coef_df_mat = as.matrix(coef_df)
powers = rep(0, nrow(coef_df))
for (row in 1:nrow(coef_df)) {
  TDC_train_features = SVM_train_features
  weights = as.vector(coef_df_mat[row, 3:(ncol(coef_df) - 1)])
  new_scores = TDC_train_features%*%weights
  SVM_labels_reordered = SVM_labels[order(-new_scores)]
  q_val = TDC_flex_c(SVM_labels_reordered == -1, SVM_labels_reordered == 1)
  power = sum(q_val <= 0.01 & SVM_labels_reordered == 1)
  powers[row] = power
}

ind.max = which.max(powers)
weights = as.vector(coef_df_mat[ind.max, 3:(ncol(coef_df) - 1)])
combined_file_features = combined_file %>% select(ExpMass, CalcMass, deltCn, XCorr, TailorScore,
                         PepLen, Charge2, Charge3, lnNumSP, dM, absdM,
)
attach(combined_file_features)
combined_file_features = combined_file_features %>% select(-Charge2, -Charge3) %>% scale() %>% as.data.frame()
combined_file_features = cbind(combined_file_features, Charge2) 
combined_file_features = as.matrix(combined_file_features)
new_scores = combined_file_features%*%weights
new_labels = combined_file$Label[order(-new_scores)]

q_val = TDC_flex_c(new_labels == -1, new_labels == 1)
power_final = sum(q_val <= 0.01 & new_labels == 1) #4356

#compare this to narrow TDC
new_labels = combined_file$Label[order(-combined_file$TailorScore)]
q_val = TDC_flex_c(new_labels == -1, new_labels == 1)
power_final = sum(q_val <= 0.01 & new_labels == 1) #4452


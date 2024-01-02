#example on running Adaptive KFs

#libraries
library(adaptiveKnockoff)
library(randomForest)
library(tidyverse)
library(gam)
source("filter_RF.R")
source("filter_EM.R")

#seed
set.seed(1012024)

#read in example data
df = read_delim('datasets/PXD019354/crux-output/open_1_0.make-pin.pin')
peptide_list = read_delim('datasets/PXD019354/index-0/tide-index.peptides.txt')

df = df[sample(nrow(df)), ]
df = df[order(df$TailorScore, decreasing = TRUE), ]
df$original_target_sequence = substring(df$Peptide, 3, nchar(df$Peptide) - 2)
df$original_target_sequence = gsub('[0-9]+|\\.|\\[|\\]', '', df$original_target_sequence)
df$original_target_sequence[df$Label == -1] = peptide_list$target[match(df$original_target_sequence[df$Label == -1], peptide_list$`decoy(s)`)]

df = df[duplicated(df$original_target_sequence), ]

W = df$TailorScore*(df$Label)
z = df %>% select('deltLCn', 'deltCn', 'Charge1', 'Charge2', 'Charge3', 'Charge4', 'Charge5',  'XCorr', 'PepLen', 'lnNumDSP', 'dM', 'absdM')
z = as.matrix(z)

props = (1:9)/10
results = rep(0, length(props))
for (i in 1:length(props)) {
  print(i)
  result = filter_RF(W, z, reveal_prop = props[i], mute = FALSE)
  results[i] = result
}

mean_time = mean(results)

print('Time for RF:')
print(paste('average time in minutes:', mean_time*nrow(df)*0.9))
print(paste('sd time in minutes:', sd(results)*nrow(df)*0.9))

W = df$TailorScore
W = W - min(W)
W = W*(df$Label)
results = rep(0, length(props))
for (i in 1:length(props)) {
  print(i)
  result = filter_EM(W, z, reveal_prop = props[i], mute = FALSE)
  results[i] = result
}

mean_time = mean(results)

print('Time for EM:')
print(paste('average time in minutes:', mean_time*nrow(df)*0.9))
print(paste('sd time in minutes:', sd(results)*nrow(df)*0.9))


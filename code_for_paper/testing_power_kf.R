#example on running Adaptive KFs

#libraries
library(adaptiveKnockoff)
library(randomForest)
library(tidyverse)
library(gam)
source("filter_RF.R")
source("filter_EM.R")
source("filter_EM_full.R")

#seed
set.seed(1012024)

#read in example data
df = read_delim('datasets/PXD022257/crux-output/open_1_0.make-pin.pin')
peptide_list = read_delim('datasets/PXD022257/index-0/tide-index.peptides.txt')

df = df[sample(nrow(df)), ]
df = df[order(df$TailorScore, decreasing = TRUE), ]
df$original_target_sequence = substring(df$Peptide, 3, nchar(df$Peptide) - 2)
df$original_target_sequence = gsub('[0-9]+|\\.|\\[|\\]', '', df$original_target_sequence)
df$original_target_sequence[df$Label == -1] = peptide_list$target[match(df$original_target_sequence[df$Label == -1], peptide_list$`decoy(s)`)]

df = df[!duplicated(df$original_target_sequence), ]

z = df %>% select('deltLCn', 'deltCn', 'Charge1', 'Charge2', 'Charge3', 'Charge4', 'Charge5',  'XCorr', 'PepLen', 'lnNumDSP', 'dM', 'absdM')
z = as.matrix(z)

W = df$TailorScore
W = W - min(W)
W = W*(df$Label)
W[abs(W) <= 1e-3] = 0

#number of Percolator-RESET discoveries
alpha = 0.01
target_peps = read.delim('datasets/PXD022257/crux-output/testing_adaptive_kf.peptides.txt')
decoy_peps = read.delim('datasets/PXD022257/crux-output/testing_adaptive_kf.decoy_peptides.txt')
min_score = min(target_peps$SVM_score[target_peps$q_val <= alpha])
targets = sum(target_peps$SVM_score >= min_score)
decoys = sum(decoy_peps$SVM_score >= min_score)
reset = targets + decoys

props = (1:9)/10
results = rep(0, length(props))
for (i in 1:length(props)) {
  print(i)
  result = filter_EM(W, z, reveal_prop = props[i], mute = FALSE)
  results[i] = result
}
relevant = floor((nrow(df) - reset)/nrow(df)*10)
mean_time = mean(results[1:relevant])

print('Time for EM:')
print(paste('average time in minutes:', mean_time*(nrow(df)*0.9 - reset) ))
print(paste('sd time in minutes:', sd(results)*(nrow(df)*0.9 - reset) ))


start = Sys.time()
results = filter_EM_full(W, z, alpha = c(0.05, 0.04, 0.03, 0.02, 0.01), reveal_prop = 0.1, mute = FALSE)
end = Sys.time()

print('Actual time:')
print(end - start)

print('Rejections:')
print(results$nrejs)

#plot values
adakf = rev(c(2888, 2763, 2668, 2462, 2264))
df_RESET = read.delim('datasets/PXD022257/crux-output/testing_adaptive_kf.peptides.txt')
reset = c(sum(df_RESET$q_val <= 0.01), sum(df_RESET$q_val <= 0.02), sum(df_RESET$q_val <= 0.03), sum(df_RESET$q_val <= 0.04), sum(df_RESET$q_val <= 0.05))



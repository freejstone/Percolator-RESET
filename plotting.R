#plots

#########################################################################################################################
library(tidyverse)

percolator_all = read.csv('results/percolator_narrow_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')


df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character(),
                    alpha = as.numeric())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    for (alpha in seq(0.01, 0.05, by = 0.01)) {
      print(i)
      super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', 'narrow_single_decoy_', i, '.peptides.txt', sep = ''))
      df_all[count, ] = c(sum(super_percolator$q_val <= alpha), i, dataset, alpha)
      count = count + 1
    }
  }
}

percolator_all = percolator_all[percolator_all$Threshold %in% seq(0.01, 0.05, by = 0.01), ]

all(percolator_all$PXID == df_all$PXID) #TRUE
all(percolator_all$Threshold == df_all$alpha) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_narrow)

df_all_summary = df_all %>% group_by(alpha, PXID) %>% summarize(median_PXID = median(ratio), q25_PXID = quantile(ratio, 0.25), q75_PXID = quantile(ratio, 0.75))

df_summary = df_all_summary %>% group_by(alpha) %>% summarize(median = median(median_PXID), q25 = quantile(median_PXID, 0.25), q75 = quantile(median_PXID, 0.75))

df_long <- df_summary %>%
  gather(key = "variable", value = "value", -alpha)

df_long$variable[df_long$variable == 'median'] = '0.5'
df_long$variable[df_long$variable == 'q25'] = '0.25'
df_long$variable[df_long$variable == 'q75'] = '0.75'
df_long$variable = factor(df_long$variable, levels = c('0.25', '0.5', '0.75'))
df_long_final = df_long
df_long_final$n_o = 'Narrow'
df_long_final$type = 'Percolator+'

ggplot() + geom_line(data = df_long, mapping = aes(x = alpha, y = value, group = variable, colour = variable)) +
  geom_abline(slope = 0, intercept = 1) + ylim(c(0.97, 1.03)) +
  labs(x = expression('FDR threshold'), y = 'Ratio of Percolator+ / Percolator \n(narrow)', colour = 'Quantile') +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) 

#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')


df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character(),
                    alpha = as.numeric())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    for (alpha in seq(0.01, 0.05, by = 0.01)) {
      print(i)
      super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', 'open_single_decoy_', i, '.peptides.txt', sep = ''))
      df_all[count, ] = c(sum(super_percolator$q_val <= alpha), i, dataset, alpha)
      count = count + 1
    }
  }
}

percolator_all = percolator_all[percolator_all$Threshold %in% seq(0.01, 0.05, by = 0.01), ]

all(percolator_all$PXID == df_all$PXID) #TRUE
all(percolator_all$Threshold == df_all$alpha) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)

df_all_summary = df_all %>% group_by(alpha, PXID) %>% summarize(median_PXID = median(ratio), q25_PXID = quantile(ratio, 0.25), q75_PXID = quantile(ratio, 0.75))

df_summary = df_all_summary %>% group_by(alpha) %>% summarize(median = median(median_PXID), q25 = quantile(median_PXID, 0.25), q75 = quantile(median_PXID, 0.75))

df_long <- df_summary %>%
  gather(key = "variable", value = "value", -alpha)

df_long$variable[df_long$variable == 'median'] = '0.5'
df_long$variable[df_long$variable == 'q25'] = '0.25'
df_long$variable[df_long$variable == 'q75'] = '0.75'
df_long$variable = factor(df_long$variable, levels = c('0.25', '0.5', '0.75'))
df_long$type = 'Percolator+'
df_long$n_o = 'Open'
df_long_final = rbind(df_long_final, df_long)

ggplot() + geom_line(data = df_long, mapping = aes(x = alpha, y = value, group = variable, colour = variable)) +
  geom_abline(slope = 0, intercept = 1) + ylim(c(0.97, 1.03)) +
  labs(x = expression('FDR threshold'), y = 'Ratio of Percolator+ / Percolator \n(open)', colour = 'Quantile') +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) 


#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_narrow_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')


df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character(),
                    alpha = as.numeric())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    for (alpha in seq(0.01, 0.05, by = 0.01)) {
      print(i)
      super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', 'narrow_extra_decoy_', i, '.peptides.txt', sep = ''))
      df_all[count, ] = c(sum(super_percolator$q_val <= alpha), i, dataset, alpha)
      count = count + 1
    }
  }
}

percolator_all = percolator_all[percolator_all$Threshold %in% seq(0.01, 0.05, by = 0.01), ]
percolator_all$power_narrow_avg = 0

for (dataset in datasets) {
  for (i in 0:9) {
    for (alpha in seq(0.01, 0.05, by = 0.01)) {
      print(i)
      percolator_all$power_narrow_avg[percolator_all$PXID == dataset & percolator_all$Threshold == alpha & percolator_all$dcy_indx == i] = (percolator_all$power_narrow[percolator_all$PXID == dataset & percolator_all$Threshold == alpha & percolator_all$dcy_indx == i] + percolator_all$power_narrow[percolator_all$PXID == dataset & percolator_all$Threshold == alpha & percolator_all$dcy_indx == ((i + 1) %% 10)])/2
      count = count + 1
    }
  }
}

all(percolator_all$PXID == df_all$PXID) #TRUE
all(percolator_all$Threshold == df_all$alpha) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_narrow_avg)

df_all_summary = df_all %>% group_by(alpha, PXID) %>% summarize(median_PXID = median(ratio), q25_PXID = quantile(ratio, 0.25), q75_PXID = quantile(ratio, 0.75))

df_summary = df_all_summary %>% group_by(alpha) %>% summarize(median = median(median_PXID), q25 = quantile(median_PXID, 0.25), q75 = quantile(median_PXID, 0.75))

df_long <- df_summary %>%
  gather(key = "variable", value = "value", -alpha)

df_long$variable[df_long$variable == 'median'] = '0.5'
df_long$variable[df_long$variable == 'q25'] = '0.25'
df_long$variable[df_long$variable == 'q75'] = '0.75'
df_long$variable = factor(df_long$variable, levels = c('0.25', '0.5', '0.75'))
df_long$type = 'Percolator extra+'
df_long$n_o = 'Narrow'
df_long_final = rbind(df_long_final, df_long)

ggplot() + geom_line(data = df_long, mapping = aes(x = alpha, y = value, group = variable, colour = variable)) +
  geom_abline(slope = 0, intercept = 1) + ylim(c(0.97, 1.03)) +
  labs(x = expression('FDR threshold'), y = 'Ratio of Percolator extra+ / Percolator \n(narrow)', colour = 'Quantile') +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) 



#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')


df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character(),
                    alpha = as.numeric())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    for (alpha in seq(0.01, 0.05, by = 0.01)) {
      print(i)
      super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', 'open_extra_decoy_', i, '.peptides.txt', sep = ''))
      df_all[count, ] = c(sum(super_percolator$q_val <= alpha), i, dataset, alpha)
      count = count + 1
    }
  }
}

percolator_all = percolator_all[percolator_all$Threshold %in% seq(0.01, 0.05, by = 0.01), ]
percolator_all$power_open_avg = 0

for (dataset in datasets) {
  for (i in 0:9) {
    for (alpha in seq(0.01, 0.05, by = 0.01)) {
      print(i)
      percolator_all$power_open_avg[percolator_all$PXID == dataset & percolator_all$Threshold == alpha & percolator_all$dcy_indx == i] = (percolator_all$power_open[percolator_all$PXID == dataset & percolator_all$Threshold == alpha & percolator_all$dcy_indx == i] + percolator_all$power_open[percolator_all$PXID == dataset & percolator_all$Threshold == alpha & percolator_all$dcy_indx == ((i + 1) %% 10)])/2
      count = count + 1
    }
  }
}

all(percolator_all$PXID == df_all$PXID) #TRUE
all(percolator_all$Threshold == df_all$alpha) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open_avg)

df_all_summary = df_all %>% group_by(alpha, PXID) %>% summarize(median_PXID = median(ratio), q25_PXID = quantile(ratio, 0.25), q75_PXID = quantile(ratio, 0.75))

df_summary = df_all_summary %>% group_by(alpha) %>% summarize(median = median(median_PXID), q25 = quantile(median_PXID, 0.25), q75 = quantile(median_PXID, 0.75))

df_long <- df_summary %>%
  gather(key = "variable", value = "value", -alpha)

df_long$variable[df_long$variable == 'median'] = '0.5'
df_long$variable[df_long$variable == 'q25'] = '0.25'
df_long$variable[df_long$variable == 'q75'] = '0.75'
df_long$variable = factor(df_long$variable, levels = c('0.25', '0.5', '0.75'))
df_long$type = 'Percolator extra+'
df_long$n_o = 'Open'
df_long_final = rbind(df_long_final, df_long)
df_long_final$type = factor(df_long_final$type, levels = c('Percolator+', 'Percolator extra+'))

ggplot() + geom_line(data = df_long_final, mapping = aes(x = alpha, y = value, group = interaction(variable, n_o, type), colour = variable, linetype = n_o)) + facet_grid(.~type) + 
  geom_abline(slope = 0, intercept = 1) + ylim(c(0.97, 1.03)) + scale_linetype_manual("Narrow/Open",values=c("Open"=2,"Narrow"=1)) +
  labs(x = expression('FDR threshold'), y = 'Ratio of Percolator extra+ / Percolator \n(open)', colour = 'Quantile') +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) 


#########################################################################################################################

datasets = list.files(path = 'datasets', pattern = 'PXD')

df_all_single = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character(),
                    alpha = as.numeric())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    for (alpha in seq(0.01, 0.05, by = 0.01)) {
      print(i)
      super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', 'narrow_single_decoy_', i, '.peptides.txt', sep = ''))
      df_all_single[count, ] = c(sum(super_percolator$q_val <= alpha), i, dataset, alpha)
      count = count + 1
    }
  }
}


df_all_extra = data.frame(power = as.numeric(),
                     dcy_index = as.numeric(),
                     PXID = as.character(),
                     alpha = as.numeric())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    for (alpha in seq(0.01, 0.05, by = 0.01)) {
      print(i)
      super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', 'narrow_extra_decoy_', i, '.peptides.txt', sep = ''))
      df_all_extra[count, ] = c(sum(super_percolator$q_val <= alpha), i, dataset, alpha)
      count = count + 1
    }
  }
}


df_all_extra_summary = df_all_extra %>% group_by(alpha, PXID) %>% summarize(variance = var(as.numeric(power)))
df_all_single_summary = df_all_single %>% group_by(alpha, PXID) %>% summarize(variance = var(as.numeric(power)))

df_all_extra_summary$ratio = df_all_extra_summary$variance/df_all_single_summary$variance

df_all_summary = df_all_extra_summary %>% group_by(alpha) %>% summarize(IQR_all = median(ratio), IQR_q25 = quantile(ratio, 0.25), IQR_q75 = quantile(ratio, 0.75))

df_long <- df_all_summary %>%
  gather(key = "variable", value = "value", -alpha)
df_long$variable[df_long$variable == 'IQR_all'] = '0.5'
df_long$variable[df_long$variable == 'IQR_q25'] = '0.25'
df_long$variable[df_long$variable == 'IQR_q75'] = '0.75'
df_long$variable = factor(df_long$variable, levels = c('0.25', '0.5', '0.75'))

ggplot() + geom_line(data = df_long, mapping = aes(x = alpha, y = log(value, 2), group = variable, colour = variable)) +
  geom_abline(slope = 0, intercept = 0) + ylim(-2, 1) +
  labs(x = expression('FDR threshold'), y = 'Ratio of Percolator extra+ variance / Percolator+ \n(narrow) (log scale)', colour = 'Quantile') +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) 

#########################################################################################################################

datasets = list.files(path = 'datasets', pattern = 'PXD')

df_all_single = data.frame(power = as.numeric(),
                           dcy_index = as.numeric(),
                           PXID = as.character(),
                           alpha = as.numeric())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    for (alpha in seq(0.01, 0.05, by = 0.01)) {
      print(i)
      super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', 'open_single_decoy_', i, '.peptides.txt', sep = ''))
      df_all_single[count, ] = c(sum(super_percolator$q_val <= alpha), i, dataset, alpha)
      count = count + 1
    }
  }
}


df_all_extra = data.frame(power = as.numeric(),
                          dcy_index = as.numeric(),
                          PXID = as.character(),
                          alpha = as.numeric())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    for (alpha in seq(0.01, 0.05, by = 0.01)) {
      print(i)
      super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', 'open_extra_decoy_', i, '.peptides.txt', sep = ''))
      df_all_extra[count, ] = c(sum(super_percolator$q_val <= alpha), i, dataset, alpha)
      count = count + 1
    }
  }
}


df_all_extra_summary = df_all_extra %>% group_by(alpha, PXID) %>% summarize(variance = var(as.numeric(power)))
df_all_single_summary = df_all_single %>% group_by(alpha, PXID) %>% summarize(variance = var(as.numeric(power)))

df_all_extra_summary$ratio = df_all_extra_summary$variance/df_all_single_summary$variance

df_all_summary = df_all_extra_summary %>% group_by(alpha) %>% summarize(IQR_all = median(ratio), IQR_q25 = quantile(ratio, 0.25), IQR_q75 = quantile(ratio, 0.75))

df_long <- df_all_summary %>%
  gather(key = "variable", value = "value", -alpha)
df_long$variable[df_long$variable == 'IQR_all'] = '0.5'
df_long$variable[df_long$variable == 'IQR_q25'] = '0.25'
df_long$variable[df_long$variable == 'IQR_q75'] = '0.75'
df_long$variable = factor(df_long$variable, levels = c('0.25', '0.5', '0.75'))

ggplot() + geom_line(data = df_long, mapping = aes(x = alpha, y = log(value, 2), group = variable, colour = variable)) +
  geom_abline(slope = 0, intercept = 0) + ylim(-2, 1) +
  labs(x = expression('FDR threshold'), y = 'Ratio of Percolator extra+ variance / Percolator+ \n(open) (log scale)', colour = 'Quantile') +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) 


#########################################################################################################################

datasets = list.files(path = 'datasets', pattern = 'PXD')

percolator_all = read.csv('results/percolator_all.csv')
names(percolator_all)[2] = "alpha"

df_all_extra = data.frame(power = as.numeric(),
                          dcy_index = as.numeric(),
                          PXID = as.character(),
                          alpha = as.numeric())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    for (alpha in seq(0.01, 0.05, by = 0.01)) {
      print(i)
      super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', 'open_extra_decoy_', i, '.peptides.txt', sep = ''))
      df_all_extra[count, ] = c(sum(super_percolator$q_val <= alpha), i, dataset, alpha)
      count = count + 1
    }
  }
}

df_all_extra_summary = df_all_extra %>% group_by(alpha, PXID) %>% summarize(variance = var(as.numeric(power)))
percolator_all = percolator_all[percolator_all$alpha %in% seq(0.01, 0.05, by = 0.01), ]
percolator_all_summary = percolator_all %>% group_by(alpha, PXID) %>% summarize(variance = var(as.numeric(power_open)))

df_all_extra_summary$ratio = df_all_extra_summary$variance/percolator_all_summary$variance

df_all_summary = df_all_extra_summary %>% group_by(alpha) %>% summarize(IQR_all = median(ratio), IQR_q25 = quantile(ratio, 0.25), IQR_q75 = quantile(ratio, 0.75))

df_long <- df_all_summary %>%
  gather(key = "variable", value = "value", -alpha)
df_long$variable[df_long$variable == 'IQR_all'] = '0.5'
df_long$variable[df_long$variable == 'IQR_q25'] = '0.25'
df_long$variable[df_long$variable == 'IQR_q75'] = '0.75'
df_long$variable = factor(df_long$variable, levels = c('0.25', '0.5', '0.75'))

ggplot() + geom_line(data = df_long, mapping = aes(x = alpha, y = value, group = variable, colour = variable)) +
  geom_abline(slope = 0, intercept = 1) + 
  labs(x = expression('FDR threshold'), y = 'Ratio of Percolator++ variance / Percolator+ \n(open)', colour = 'Quantile') +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10)) 

  
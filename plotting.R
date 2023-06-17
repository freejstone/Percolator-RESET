#plots

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')

df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  for (i in 0:9) {
    super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', i, '.power.txt', sep = ''))
    df_all[count, ] = c(super_percolator$true_power[5], i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_summary1 = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))

#########################################################################################################################

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')

df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  for (i in 0:9) {
    super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', i, '_top_positive.power.txt', sep = ''))
    df_all[count, ] = c(super_percolator$true_power[5], i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_summary = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))

#########################################################################################################################


percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')

datasets = setdiff(datasets, c('PXD019354'))

df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', 'extra_decoy_', i, '_SVM.peptides.txt', sep = ''))
    df_all[count, ] = c(sum(super_percolator$q_vals[super_percolator$Label == 1] <= 0.01), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]

percolator_all = percolator_all[percolator_all$PXID != 'PXD019354', ]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_summary = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))




#########################################################################################################################


percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')


df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = read.delim(paste('datasets/', dataset, '/crux-output/', 'adaptive_', i, '_grid30.power.txt', sep = ''))
    df_all[count, ] = c(super_percolator$true_power[10], i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_summary = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))


#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')


df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'stratified_', i, '.log.txt', sep = ''))
    indx = which(grepl('iteration: 4', super_percolator))
    indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_summary_strat = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))



#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')


df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'no_stratified_', i, '.log.txt', sep = ''))
    indx = which(grepl('iteration: 4', super_percolator))
    indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_summary_no_strat_50 = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))



#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')


df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'no_stratified_pinit_025_', i, '.log.txt', sep = ''))
    indx = which(grepl('iteration: 4', super_percolator))
    indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_summary_no_strat_25 = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))



#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')


df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'no_stratified_pinit_075_', i, '.log.txt', sep = ''))
    indx = which(grepl('iteration: 4', super_percolator))
    indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_summary_no_strat_75 = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))



#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')


df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'stratified_freq_', i, '.log.txt', sep = ''))
    indx = which(grepl('iteration: 4', super_percolator))
    indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_summary_strat_freq = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))

#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')
datasets = datasets[datasets != "PXD019354"]

df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'extra_decoy_', i, '_PSM_level.log.txt', sep = ''))
    indx = which(grepl('peptides discovered', super_percolator))
    #indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]
percolator_all = percolator_all[percolator_all$PXID != "PXD019354", ]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_summary_extra_decoy_psm = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))

#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')
#datasets = datasets[datasets != "PXD019354"]

df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'extra_decoy_', i, '_pep_level_075.log.txt', sep = ''))
    indx = which(grepl('iteration: 4', super_percolator))
    indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]
#percolator_all = percolator_all[percolator_all$PXID != "PXD019354", ]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_summary_extra_decoy_pep = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))

#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_narrow_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')
#datasets = datasets[datasets != "PXD019354"]

df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'no_stratified_narrow_', i, '.log.txt', sep = ''))
    indx = which(grepl('iteration: 4', super_percolator))
    indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]
#percolator_all = percolator_all[percolator_all$PXID != "PXD019354", ]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_narrow)
df_all_summary_narrow = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))


#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')
#datasets = datasets[datasets != "PXD019354"]

df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'extra_decoy_', i, '_pep_level_075.log.txt', sep = ''))
    indx = which(grepl('iteration: 4', super_percolator))
    indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]
#percolator_all = percolator_all[percolator_all$PXID != "PXD019354", ]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_extra_pep_075 = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))


#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')
#datasets = datasets[datasets != "PXD019354"]

df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'extra_decoy_', i, '_pep_level_050.log.txt', sep = ''))
    indx = which(grepl('iteration: 4', super_percolator))
    indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]
#percolator_all = percolator_all[percolator_all$PXID != "PXD019354", ]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_extra_pep_050 = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))


#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_narrow_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')
#datasets = datasets[datasets != "PXD019354"]

df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'extra_decoy_', i, '_pep_level_narrow_075.log.txt', sep = ''))
    indx = which(grepl('iteration: 4', super_percolator))
    indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]
#percolator_all = percolator_all[percolator_all$PXID != "PXD019354", ]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_narrow)
df_all_extra_pep_075_narrow = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))



#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_narrow_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')
datasets = datasets[datasets != "PXD019354"]

df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'extra_decoy_', i, '_PSM_level_narrow.log.txt', sep = ''))
    indx = which(grepl('peptides discovered', super_percolator))
    #indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]
percolator_all = percolator_all[percolator_all$PXID != "PXD019354", ]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_narrow)
df_all_extra_psm_narrow = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))


#########################################################################################################################
library(tidyr)

percolator_all = read.csv('results/percolator_all.csv')

datasets = list.files(path = 'datasets', pattern = 'PXD')
datasets = datasets[datasets != "PXD019354"]

df_all = data.frame(power = as.numeric(),
                    dcy_index = as.numeric(),
                    PXID = as.character())
count = 1
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    print(i)
    super_percolator = readLines(paste('datasets/', dataset, '/crux-output/', 'extra_decoy_', i, '_PSM_level.log.txt', sep = ''))
    indx = which(grepl('peptides discovered', super_percolator))
    #indx = indx[length(indx)] + 1
    df_all[count, ] = c(parse_number(super_percolator[indx]), i, dataset)
    count = count + 1
  }
}

percolator_all = percolator_all[percolator_all$Threshold == 0.01,]
percolator_all = percolator_all[percolator_all$PXID != "PXD019354", ]

all(percolator_all$PXID == df_all$PXID) #TRUE

df_all$ratio = as.numeric(df_all$power)/as.numeric(percolator_all$power_open)
df_all_extra_psm_open = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))

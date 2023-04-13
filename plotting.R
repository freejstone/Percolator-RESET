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
df_all_summary = df_all %>% group_by(PXID) %>% summarize(median = median(ratio), q25 = quantile(ratio, 0.25), q75 = quantile(ratio, 0.75))
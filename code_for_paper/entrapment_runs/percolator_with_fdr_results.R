#percolator with fdr entrapments

#parameters
sizes = c('full', 'half', '075', '025')
alpha_vec = seq(0.01, 0.1, by = 0.01)
d_vec = c(0:199)
count = 1

df_all = data.frame(Threshold = numeric(),
                    power = numeric(),
                    est_FDP = numeric(),
                    index = numeric(),
                    sample_size = numeric()
                    )

for (size in sizes) {
  for (d in d_vec) {
    print(d)
    CONGA_df = read.delim(paste('percolator_with_fdr_results/', 'open_single_decoy_', size, '_', d, '_ind_fix.peptides.txt', sep = ''))
    CONGA_df = CONGA_df[CONGA_df$Label == 1, ]
    for (alpha in alpha_vec) {
      power = sum(CONGA_df$q_val <= alpha)
      if (power > 0) {
        est_FDP = sum(grepl('ENTR', CONGA_df$Proteins[CONGA_df$q_val <= alpha]))/power
      } else {
        est_FDP = 0
      }
      df_all[count,] = c(alpha,  power, est_FDP, d, size)
      count = count + 1
    }
  }
}

write.csv(df_all, 'results/percolator_with_fdr_open_entrapments_ind_fix.csv')
####################################################################################

#percolator with fdr entrapments

#parameters
sizes = c('full', 'half', '075', '025')
alpha_vec = seq(0.01, 0.1, by = 0.01)
d_vec = c(0:199)
count = 1

df_all = data.frame(Threshold = numeric(),
                    power = numeric(),
                    est_FDP = numeric(),
                    index = numeric(),
                    sample_size = numeric()
)

for (size in sizes) {
  for (d in d_vec) {
    print(d)
    CONGA_df = read.delim(paste('percolator_with_fdr_results/', 'narrow_single_decoy_', size, '_', d, '_ind_fix.peptides.txt', sep = ''))
    CONGA_df = CONGA_df[CONGA_df$Label == 1, ]
    for (alpha in alpha_vec) {
      power = sum(CONGA_df$q_val <= alpha)
      if (power > 0) {
        est_FDP = sum(grepl('ENTR', CONGA_df$Proteins[CONGA_df$q_val <= alpha]))/power
      } else {
        est_FDP = 0
      }
      df_all[count,] = c(alpha,  power, est_FDP, d, size)
      count = count + 1
    }
  }
}

write.csv(df_all, 'results/percolator_with_fdr_narrow_entrapments_ind_fix.csv')
####################################################################################




#percolator with fdr entrapments

#parameters
sizes = c('full', 'half', '075', '025')
alpha_vec = seq(0.01, 0.1, by = 0.01)
d_vec = c(0:199)
count = 1

df_all = data.frame(Threshold = numeric(),
                    power = numeric(),
                    est_FDP = numeric(),
                    index = numeric(),
                    sample_size = numeric()
)

for (size in sizes) {
  for (d in d_vec) {
    print(d)
    CONGA_df = read.delim(paste('percolator_with_fdr_results/', 'open_extra_decoy_', size, '_', d, '_ind_fix.peptides.txt', sep = ''))
    CONGA_df = CONGA_df[CONGA_df$Label == 1, ]
    for (alpha in alpha_vec) {
      power = sum(CONGA_df$q_val <= alpha)
      if (power > 0) {
        est_FDP = sum(grepl('ENTR', CONGA_df$Proteins[CONGA_df$q_val <= alpha]))/power
      } else {
        est_FDP = 0
      }
      df_all[count,] = c(alpha,  power, est_FDP, d, size)
      count = count + 1
    }
  }
}

write.csv(df_all, 'results/percolator_with_fdr_open_extra_entrapments_ind_fix.csv')
####################################################################################

#percolator with fdr entrapments

#parameters
sizes = c('full', '075', 'half', '025')
alpha_vec = seq(0.01, 0.1, by = 0.01)
d_vec = c(0:199)
count = 1

df_all = data.frame(Threshold = numeric(),
                    power = numeric(),
                    est_FDP = numeric(),
                    index = numeric(),
                    sample_size = numeric()
)

for (size in sizes) {
  for (d in d_vec) {
    print(d)
    CONGA_df = read.delim(paste('percolator_with_fdr_results/', 'narrow_extra_decoy_', size, '_', d, '_ind_fix.peptides.txt', sep = ''))
    CONGA_df = CONGA_df[CONGA_df$Label == 1, ]
    for (alpha in alpha_vec) {
      power = sum(CONGA_df$q_val <= alpha)
      if (power > 0) {
        est_FDP = sum(grepl('ENTR', CONGA_df$Proteins[CONGA_df$q_val <= alpha]))/power
      } else {
        est_FDP = 0
      }
      df_all[count,] = c(alpha,  power, est_FDP, d, size)
      count = count + 1
    }
  }
}

write.csv(df_all, 'results/percolator_with_fdr_narrow_extra_entrapments_ind_fix.csv')
####################################################################################



#parameters
sizes = c('full', 'half', '075', '025')
alpha_vec = seq(0.01, 0.1, by = 0.01)
d_vec = c(0:199)
count = 1

df_all = data.frame(Threshold = numeric(),
                    power = numeric(),
                    est_FDP = numeric(),
                    index = numeric(),
                    sample_size = numeric()
)

for (size in sizes) {
  for (d in d_vec) {
    print(d)
    CONGA_df = read.delim(paste('percolator_with_fdr_results/', 'open_single_decoy_', size, '_', d, '_no_pair_fix.peptides.txt', sep = ''))
    CONGA_df = CONGA_df[CONGA_df$Label == 1, ]
    for (alpha in alpha_vec) {
      power = sum(CONGA_df$q_val <= alpha)
      if (power > 0) {
        est_FDP = sum(grepl('ENTR', CONGA_df$Proteins[CONGA_df$q_val <= alpha]))/power
      } else {
        est_FDP = 0
      }
      df_all[count,] = c(alpha,  power, est_FDP, d, size)
      count = count + 1
    }
  }
}

write.csv(df_all, 'results/percolator_with_fdr_open_entrapments_no_pair_fix.csv')
####################################################################################

#percolator with fdr entrapments

#parameters
sizes = c('full', 'half', '075', '025')
alpha_vec = seq(0.01, 0.1, by = 0.01)
d_vec = c(0:199)
count = 1

df_all = data.frame(Threshold = numeric(),
                    power = numeric(),
                    est_FDP = numeric(),
                    index = numeric(),
                    sample_size = numeric()
)

for (size in sizes) {
  for (d in d_vec) {
    print(d)
    CONGA_df = read.delim(paste('percolator_with_fdr_results/', 'narrow_single_decoy_', size, '_', d, '_no_pair_fix.peptides.txt', sep = ''))
    CONGA_df = CONGA_df[CONGA_df$Label == 1, ]
    for (alpha in alpha_vec) {
      power = sum(CONGA_df$q_val <= alpha)
      if (power > 0) {
        est_FDP = sum(grepl('ENTR', CONGA_df$Proteins[CONGA_df$q_val <= alpha]))/power
      } else {
        est_FDP = 0
      }
      df_all[count,] = c(alpha,  power, est_FDP, d, size)
      count = count + 1
    }
  }
}

write.csv(df_all, 'results/percolator_with_fdr_narrow_entrapments_no_pair_fix.csv')
####################################################################################




#percolator with fdr entrapments

#parameters
sizes = c('full', 'half', '075', '025')
alpha_vec = seq(0.01, 0.1, by = 0.01)
d_vec = c(0:199)
count = 1

df_all = data.frame(Threshold = numeric(),
                    power = numeric(),
                    est_FDP = numeric(),
                    index = numeric(),
                    sample_size = numeric()
)

for (size in sizes) {
  for (d in d_vec) {
    print(d)
    CONGA_df = read.delim(paste('percolator_with_fdr_results/', 'open_extra_decoy_', size, '_', d, '_no_pair_fix.peptides.txt', sep = ''))
    CONGA_df = CONGA_df[CONGA_df$Label == 1, ]
    for (alpha in alpha_vec) {
      power = sum(CONGA_df$q_val <= alpha)
      if (power > 0) {
        est_FDP = sum(grepl('ENTR', CONGA_df$Proteins[CONGA_df$q_val <= alpha]))/power
      } else {
        est_FDP = 0
      }
      df_all[count,] = c(alpha,  power, est_FDP, d, size)
      count = count + 1
    }
  }
}

write.csv(df_all, 'results/percolator_with_fdr_open_extra_entrapments_no_pair_fix.csv')
####################################################################################

#percolator with fdr entrapments

#parameters
sizes = c('full', '075', 'half', '025')
alpha_vec = seq(0.01, 0.1, by = 0.01)
d_vec = c(0:199)
count = 1

df_all = data.frame(Threshold = numeric(),
                    power = numeric(),
                    est_FDP = numeric(),
                    index = numeric(),
                    sample_size = numeric()
)

for (size in sizes) {
  for (d in d_vec) {
    print(d)
    CONGA_df = read.delim(paste('percolator_with_fdr_results/', 'narrow_extra_decoy_', size, '_', d, '_no_pair_fix.peptides.txt', sep = ''))
    CONGA_df = CONGA_df[CONGA_df$Label == 1, ]
    for (alpha in alpha_vec) {
      power = sum(CONGA_df$q_val <= alpha)
      if (power > 0) {
        est_FDP = sum(grepl('ENTR', CONGA_df$Proteins[CONGA_df$q_val <= alpha]))/power
      } else {
        est_FDP = 0
      }
      df_all[count,] = c(alpha,  power, est_FDP, d, size)
      count = count + 1
    }
  }
}

write.csv(df_all, 'results/percolator_with_fdr_narrow_extra_entrapments_no_pair_fix.csv')
####################################################################################



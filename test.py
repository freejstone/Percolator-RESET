#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  7 20:59:52 2023

@author: jfre0619
"""

p = 1/4

#get second sample proportion
p_next = 1 - 1/(2*(1 - p))

#do scale
df_all = spf.do_scale(df_all)

q_val = uf.TDC_flex_c(
    df_all.Label == -1, df_all.Label == 1)

df_nah_temp = df_all[~(q_val <= 0.3)].copy()

df_new_temp = df_all[(q_val <= 0.3)].copy()



#create target-decoys at pseudolevel
train_decoys1 = df_new_temp.loc[(df_new_temp['Label']
                           == -1)].sample(frac=p).copy()

train_decoys2 = df_nah_temp

train_decoys = pd.concat([train_decoys1, train_decoys2])
#train_decoys = df_all[df_all['Label'] == -1].iloc[1::2, :].copy()


 #create train dataframe
#train_decoy_indxs = random.choices(
#    [True, False], k=sum(df_all['Label'] == -1), weights = [p, 1-p])
#train_decoys = df_all[df_all['Label'] == -1].copy()
#train_decoys = train_decoys[train_decoy_indxs]

train_power, std_power, true_power, df_new, train_decoys_new = do_svm(df_new_temp, train_decoys, folds=3, Cs=[
            0.1, 1, 10], p = p, total_iter=10, kernel=kernel, alpha=FDR_threshold, train_alpha=0.01, degree=2, remove=remove, top_positive=top_positive)

#train_decoys_new, df_new = do_lda(df_new_temp, train_decoys, total_iter=1, p=p, alpha=0.01, train_alpha=0.01, remove=remove, top_positive=True)

#go to to the top 20% and repeat
df_new['trained'] = 0
train_decoys_new['trained'] = 1
df_new_temp = pd.concat([df_new, train_decoys_new])

#do scaling again
df_new_temp = spf.do_scale(df_new_temp)

df_new_temp = df_new_temp.sort_values(
    by='SVM_score', ascending=False).reset_index(drop=True)

q_val = uf.TDC_flex_c(
    df_new_temp.Label == -1, df_new_temp.Label == 1)

df_nah_temp = df_new_temp[~(q_val <= 0.1)].copy()

df_new_temp = df_new_temp[(q_val <= 0.1)]



train_decoys1 = df_new_temp[df_new_temp.trained > 0 ]

train_decoys2 = df_new_temp[(df_new_temp['Label']
                           == -1) & (df_new_temp.trained < 0)].sample(frac=p_next).copy()

train_decoys3 = df_nah_temp

df_new_temp = df_new_temp.drop(['trained'], axis = 1)

train_decoys_new = pd.concat([train_decoys1, train_decoys2, train_decoys3])

train_decoys_new = train_decoys_new.drop(['trained'], axis = 1)

#do_lda(df_new_temp, train_decoys, total_iter=10, p=1/2, alpha=0.01, train_alpha=0.01, remove=remove, top_positive=True)


asdf1, asdf2, asdf3, asdf4, asdf5 = do_svm(df_new_temp, train_decoys_new, folds=folds, Cs=[
            0.1, 1, 10], p =1/2, total_iter=10, kernel=kernel, alpha=FDR_threshold, train_alpha=0.01, degree=degree, remove=remove, top_positive=top_positive)




train_power.append(train_power_next)
std_power.append(std_power_next)
true_power.append(true_power_next)
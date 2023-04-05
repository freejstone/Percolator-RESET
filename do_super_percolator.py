#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 17:15:01 2023

@author: jackfreestone

This module performs the percolator algorithm with strict FDR control 
"""

import numpy as np
import pandas as pd
from sklearn import svm
from sklearn.model_selection import cross_validate
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import GridSearchCV

#########################################################################################################
def TDC_flex_c(decoy_wins, target_wins, BC1 = 1, c = 1/2, lam = 1/2):
    '''
    Parameters
    ----------
    decoy_wins : boolean series
        indicates the position of the decoy wins.
    target_wins : booolean series
        indicates the position of the target wins.
    BC1 : float, optional
        The penalty applied to the estimated number of decoys. The default is 1.
    c : float, optional
        probability of a target win. The default is 1/2.
    lam : float, optional
        1 - probability of a decoy win. The default is 1/2.

    Returns
    -------
    q-values.

    '''
    nTD = np.cumsum(target_wins)
    nDD = np.cumsum(decoy_wins)
    fdps = np.minimum(1, ((BC1 + nDD)/ nTD) * (c / (1-lam)))
    qvals = fdps[::-1].cummin()[::-1]
    return(qvals)
#########################################################################################################
def peptide_level(narrow_file, open_file, peptide_list, score = 'TailorScore', precursor_bin_width = 1.0005079/4):
    '''
    Parameters
    ----------
    narrow_file : string
        path to filtered narrow pin file.
    open_file : string
        path to filtered open pin file.
    peptide_list : string
        path to target decoy pairing.
    score : string, optional
        Score function for initial direction. The default is 'TailorScore'.

    Returns
    -------
    .

    '''
    narrow_df = pd.read_table(narrow_file) #reading
    open_df = pd.read_table(open_file) #reading
    narrow_df['n_o'] = 1 #giving narrow-open column
    open_df['n_o'] = 0
    
    df_all = narrow_df.merge(open_df) #combine
    
    
    df_all['rank'] = df_all['SpecId'].apply(lambda x: x[-1]) #getting the rank
    df_all = df_all[df_all.rank <= 2]
    df_all['rank'][df_all.rank == 1] = 0
    df_all['rank'][df_all.rank == 2] = 1
    
    
    df_all = df_all.sample(frac = 1) #break ties randomly 
    df_all = df_all.sort_values(by = score, ascending = False) #sort by score
    
    peptide_list_df = pd.read_table(peptide_list) #reading
    
    df_all['original_target'] = df_all['Peptide'] #getting pairing
    df_all['original_target'] = df_all['original_target'].apply(lambda x: x[2:(len(x) - 2)])
    df_all_sub = df_all[df_all.Label == -1].copy()
    peptide_list_df.rename(columns = {'target': 'original_target', 'decoy':'Peptide'}, inplace = True)
    df_all_sub = df_all_sub.merge(peptide_list_df[['original_target', 'Peptide']], how = 'left', on = 'Peptide')
    df_all.loc[df_all.Label == -1, 'original_target'] = df_all_sub['original_target'].tolist()
    df_all = df_all.drop_duplicates(subset = 'original_target')
    df_all.pop('original_target')
    df_all.pop('enzInt')
    
    delta_mass_max = max(abs(df_all.ExpMass - df_all.CalcMass)) #binning delta masses
    breaks_p = np.arange(0, delta_mass_max + 2*precursor_bin_width, precursor_bin_width) - precursor_bin_width/2
    breaks_n = list(reversed(-breaks_p))
    breaks = pd.Series(breaks_n[0:(len(breaks_n) - 1)] + list(breaks_p), name = 'bins')
    digitized = np.digitize(df_all.ExpMass - df_all.CalcMass, breaks)
    df_all['bins'] = digitized
   
    bin_freq = df_all['bins'].value_counts() #getting bin frequencies
    bin_freq.columns = ['bin', 'freq']
    df_all = df_all.merge(bin_freq, how = 'left')
    
    return(df_all)
#########################################################################################################
def train_cv(labels, df, folds = 3, Cs = [0.1, 1, 10], kernel = 'linear', degree = 2, alpha = 0.01):
    '''
    Parameters
    ----------
    labels : panda series
        indicating the target-decoy wins.
    df : panda dataframe
        features to be trained on.
    folds : integer
        number of folds for k-fold cross-validation.
    Cs : list, optional
        the class_weights for SVM. The default is [0.1, 1, 10].
    kernel : string, optional
        the type of SVM kernel. The default is 'linear'.
    degree : integer, optional
        the degree used if polynomial is specified. The default is NULL.
    alpha : float, optional
        FDR threshold. The default is 0.01.

    Returns
    -------
    The optimal choice of parameters based off k-fold cross-validation.

    '''
    
    class_weight = [{-1: C_neg, 1: C_pos} for C_neg in Cs for C_pos in Cs if C_neg >= C_pos]
    gamma_range = np.logspace(-9, 3, 10)
    coef0_range = np.logspace(-9, 3, 10)
    
    param_grid = dict(gamma=gamma_range, coef0=coef0_range, class_weight = class_weight)
    
    cv = StratifiedShuffleSplit(n_splits=3, random_state=0)

    grid = GridSearchCV(svm.SVC(), param_grid=param_grid, cv=cv, scoring = 'add scoring here')


    
    for j in range(len(Cs)):
        for i in range(j, len(Cs)):
            C_neg, C_pos = Cs[i], Cs[j]
            
            
            if kernel in ['linear', 'sigmoid']:
                coef0s = [0, 0.5, 1]
                gamma = [1 / (df.shape[1]]
    
            clf = svm.SVC(kernel=kernel, class_weight = class_weight)
    
    
    return(results)    

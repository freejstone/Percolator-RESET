#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 17:15:01 2023

@author: jackfreestone

This module performs the percolator algorithm with strict FDR control 
"""

import numpy as np
import pandas as pd
import random
import utility_functions as uf
from sklearn import svm
from sklearn.metrics import make_scorer
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler

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
    pandas dataframe of the competed peptides.

    '''
    narrow_df = uf.read_pin(narrow_file) #reading
    open_df = uf.read_pin(open_file) #reading
    narrow_df['n_o'] = 1 #giving narrow-open column
    open_df['n_o'] = 0
    
    df_all = pd.concat([narrow_df, open_df]) #combine
    
    
    df_all['rank'] = df_all['SpecId'].apply(lambda x: int(x[-1])) #getting the rank
    df_all = df_all[df_all['rank'] <= 2]
    df_all.loc[df_all['rank'] == 1, 'rank'] = 0
    df_all.loc[df_all['rank'] == 2, 'rank'] = 1
    
    
    df_all = df_all.sample(frac = 1) #break ties randomly 
    df_all = df_all.sort_values(by = score, ascending = False) #sort by score
    
    peptide_list_df = pd.read_table(peptide_list) #reading
    
    df_all['Peptide'] = df_all['Peptide'].apply(lambda x: x[2:(len(x) - 2)])
    df_all_sub = df_all[df_all.Label == -1].copy()
    peptide_list_df.rename(columns = {'target': 'original_target', 'decoy':'Peptide'}, inplace = True)
    df_all_sub = df_all_sub.merge(peptide_list_df[['original_target', 'Peptide']], how = 'left', on = 'Peptide')
    df_all['original_target'] = df_all['Peptide']
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
    bin_freq = bin_freq.reset_index()
    bin_freq.columns = ['bins', 'freq']
    df_all = df_all.merge(bin_freq, how = 'left', on = 'bins')
    df_all.pop('bins')
    
    return(df_all)
#########################################################################################################
def custom_accuracy(y, y_pred, alpha = 0.01):
    '''
    Parameters
    ----------
    y : float
        the unobserved target-decoy labels.
    y_pred : float
        the new scores predicted for the unobserved test labels.

    Returns
    -------
    the number of discoveries in y sorted according to y_pred.

    '''
    y = [x for _, x in sorted(zip(y_pred, y), reverse = True)] #sorting y according to the predicted scores
    decoy_wins = [i == -1 for i in y]
    target_wins = [i == 1 for i in y]
    qvals = uf.TDC_flex_c(decoy_wins, target_wins, BC1 = 1, c = 3/4, lam = 3/4) #doing TDC
    disc = sum((qvals[i] <= alpha) and (y[i] == 1) for i in range(len(y))) #getting discoveries
    return(disc)

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
    the optimal choice of parameters based off k-fold cross-validation.

    '''
    
    if not kernel in ['linear', 'poly']:
        raise ValueError("kernel type not accepted.")
    
    class_weight = [{-1: C_neg, 1: C_pos} for C_neg in Cs for C_pos in Cs if C_neg >= C_pos]
    
    if kernel == 'linear':
        param_grid = dict(class_weight = class_weight)
    if kernel == 'poly':
        gamma_range = np.logspace(-9, 3, 10)
        coef0_range = np.logspace(-9, 3, 10)    
        param_grid = dict(gamma=gamma_range, coef0=coef0_range, class_weight = class_weight)
    
    my_scorer = make_scorer(custom_accuracy, alpha = alpha, greater_is_better=True, needs_threshold = True)
    
    grid = GridSearchCV(svm.SVC(), param_grid=param_grid, cv=folds, scoring = my_scorer)

    grid.fit(df, labels)
    
    pos = np.argmax(grid.cv_results_['mean_test_score'])
    
    selected_params = grid.cv_results_['params'][pos]
    
    return(selected_params)    
#########################################################################################################
def do_iterative_svm_cv(df_all, folds = 3, score = 'TailorScore', Cs = c(0.1, 1, 10), total_iter = 10, kernel = 'linear', alpha = 0.01, train_alpha = 0.01, degree = None)
    #create train dataframe
    train_decoy_indxs = random.choices([True, False], k = sum(df_all['Label'] == -1))
    train_decoys = df_all[df_all['Label'] == -1].copy()
    train_decoys = train_decoys[train_decoy_indxs]
    
    train_targets = df_all[~(df_all.index.isin(train_decoys.index))].copy()
    train_targets.loc[:, 'Label'] = 1
    
    train_df = pd.concat([train_decoys, train_targets]).reset_index(drop = True)
    train_df.sample(frac = 1)
    train_df = train_df.sort_values(by = score, ascending = False)
    
    #real df
    real_df = df_all[~(df_all.index.isin(train_decoys.index))].copy()
    
    #Preprocess dataframe
    SVM_train_labels = train_df['Label'].copy()
    SVM_train_features = train_df.copy()
    SVM_train_features.drop(['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins'],axis = 1, inplace = True)
    sds = SVM_train_features.apply(np.std, axis = 0)
    SVM_train_features = SVM_train_features.loc[:, abs(sds) > 1e-10]
    
    #scale non-binary features
    scale = StandardScaler()
    SVM_train_features.loc[:,~(SVM_train_features.columns.isin(["Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzN", "enzC", "rank", "n_o"]))] = scale.fit_transform(SVM_train_features.loc[:,~(SVM_train_features.columns.isin(["Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzN", "enzC", "rank", "n_o"]))])
    
    #getting initial positive and negative set
    q_vals = TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1, c = 3/4, lam = 3/4)
    positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= alpha)
    negative_set_indxs = (SVM_train_labels == -1)
        
    SVM_train_features_iter = SVM_train_features[positive_set_indxs | negative_set_indxs, ]
    SVM_train_labels_iter = SVM_train_labels[positive_set_indxs | negative_set_indxs]
    
    
df = df_all.copy()
labels = df_all.Label.copy()
df.drop(['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins'], axis = 1, inplace = True)






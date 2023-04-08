#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 23:47:06 2023

@author: jackfreestone

super_percolator main functions
"""
import numpy as np
import pandas as pd
import random
import utility_functions as uf
import logging
import sys
from sklearn import svm
from sklearn.metrics import make_scorer
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
#########################################################################################################
def peptide_level(df_all, peptide_list_df, precursor_bin_width = 1.0005079/4):
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
    
    #select only top 1 narrows or top 2 open psms
    df_all['rank'] = df_all['rank'].astype(int)
    df_all = df_all[((df_all['rank'] <= 2) & (df_all['n_o'] == 0)) | ((df_all['rank'] == 1) & (df_all['n_o'] == 1))]
    df_all.loc[df_all['rank'] == 1, 'rank'] = 0
    df_all.loc[df_all['rank'] == 2, 'rank'] = 1    
    
    df_all = df_all.sample(frac = 1) #break ties randomly 
    
    if 'TailorScore' in df_all.columns:
        df_all = df_all.sort_values(by = 'TailorScore', ascending = False).reset_index(drop = True) #sort by score
    else:
        df_all = df_all.sort_values(by = 'XCorr', ascending = False).reset_index(drop = True) #sort by score
    
    df_all = df_all.drop_duplicates(subset = 'Peptide') #getting best score for each Peptide
    
    df_all_sub = df_all[df_all.Label == -1].copy()
    peptide_list_df.rename(columns = {'target': 'original_target', 'decoy':'Peptide'}, inplace = True)
    df_all_sub = df_all_sub.merge(peptide_list_df[['original_target', 'Peptide']], how = 'left', on = 'Peptide')
    df_all['original_target'] = df_all['Peptide']
    df_all.loc[df_all.Label == -1, 'original_target'] = df_all_sub['original_target'].tolist()
    
    #adding both target-decoy scores so that they can be trained on as well
    if 'TailorScore' in df_all.columns:
        df_all = df_all.assign(min_tailor_score = df_all.groupby('original_target').TailorScore.transform(lambda x: min(x) if min(x) != max(x) else 0))
    if 'XCorr' in df_all.columns:
        df_all = df_all.assign(min_xcorr_score = df_all.groupby('original_target').XCorr.transform(lambda x: min(x) if min(x) != max(x) else 0))
        
    df_all = df_all.drop_duplicates(subset = 'original_target')
    df_all.drop(['original_target', 'enzInt'], axis = 1, inplace = True)
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
def train_cv(labels, df, folds = 3, Cs = [0.1, 1, 10], kernel = 'rbf', degree = 2, alpha = 0.01):
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
    
    if not kernel in ['linear', 'poly', 'rbf']:
        raise ValueError("kernel type not accepted.")
    
    class_weight = [{-1: C_neg, 1: C_pos} for C_neg in Cs for C_pos in Cs if C_neg >= C_pos]
    
    if kernel == 'linear':
        param_grid = dict(class_weight = class_weight)
    if kernel == 'rbf':
        gamma_range = [0.1, 1, 10]
        param_grid = dict(gamma=gamma_range, class_weight = class_weight)
    if kernel == 'poly':
        gamma_range = [0.1, 1, 10]
        coef0_range = [1]
        param_grid = dict(gamma=gamma_range, coef0 = coef0_range, class_weight = class_weight)
        
    
    my_scorer = make_scorer(custom_accuracy, alpha = alpha, greater_is_better=True, needs_threshold = True)
    
    grid = GridSearchCV(svm.SVC(), param_grid=param_grid, cv=folds, scoring = my_scorer)

    grid.fit(df, labels)
    
    return(grid)    
#########################################################################################################
def do_iterative_svm_cv(df_all, folds = 3, Cs = [0.1, 1, 10], total_iter = 10, kernel = 'linear', alpha = 0.01, train_alpha = 0.01, degree = None):
    df_all.drop(['min_tailor_score', 'min_xcorr_score'], axis = 1, inplace = True)
    
    #create train dataframe
    train_decoy_indxs = random.choices([True, False], k = sum(df_all['Label'] == -1))
    train_decoys = df_all[df_all['Label'] == -1].copy()
    train_decoys = train_decoys[train_decoy_indxs]
    
    train_targets = df_all[~(df_all.index.isin(train_decoys.index))].copy()
    train_targets.loc[:, 'Label'] = 1
    
    train_df = pd.concat([train_decoys, train_targets]).reset_index(drop = True)
    train_df = train_df.sample(frac = 1).reset_index(drop = True)
    pi_0s = train_df.groupby('bins').apply(lambda x: sum(x.Label == 1)/x.freq)
    pi_0s = pi_0s.reset_index()
    pi_0s.drop('level_1', axis = 1, inplace = True)
    pi_0s.columns = ['bins', 'pi_0']
    pi_0s.drop_duplicates(inplace = True)
    train_df = train_df.merge(pi_0s, how = 'left', on = 'bins')
    train_df.drop(['freq', 'bins'], axis = 1, inplace = True)
    if 'TailorScore' in train_df.columns:
        train_df = train_df.sort_values(by = 'TailorScore', ascending = False).reset_index(drop = True)
    else:
        train_df = train_df.sort_values(by = 'XCorr', ascending = False).reset_index(drop = True)
    
    #real df
    real_df = df_all[~(df_all.index.isin(train_decoys.index))].copy().reset_index(drop = True)
    real_df = real_df.merge(pi_0s, how = 'left', on = 'bins')
    real_df.drop(['freq', 'bins'], axis = 1, inplace = True)
    
    
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
    q_vals = uf.TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1, c = 3/4, lam = 3/4)
    positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= alpha)
    negative_set_indxs = (SVM_train_labels == -1)
        
    SVM_train_features_iter = SVM_train_features.loc[positive_set_indxs | negative_set_indxs, :].copy()
    SVM_train_labels_iter = SVM_train_labels.loc[positive_set_indxs | negative_set_indxs].copy()
    
    train_power, train_std, true_power = [0]*total_iter, [0]*total_iter, [0]*total_iter
    
    logging.info("Conducting iterative SVM.")
    sys.stderr.write("Conducting iterative SVM.\n")
    
    for iterate in range(total_iter):
        logging.info("iteration: %s." %(iterate))
        sys.stderr.write("iteration: %s.\n" %(iterate))
        #determining best direction with cross validation for parameter selection
        grid = train_cv(SVM_train_labels_iter, SVM_train_features_iter, folds = folds, Cs = Cs, kernel = kernel, degree = degree, alpha = train_alpha)
        best_train_power = max(grid.cv_results_['mean_test_score'])
        best_train_std = max(grid.cv_results_['std_test_score'])
        train_power[iterate] = best_train_power
        train_std[iterate] = best_train_std
        
        #the new direction
        new_scores = grid.decision_function(SVM_train_features)
        
        new_idx = pd.Series(new_scores).sort_values(ascending = False).index
        
        SVM_train_features = SVM_train_features.loc[new_idx].reset_index(drop = True)
        
        SVM_train_labels = SVM_train_labels.loc[new_idx].reset_index(drop = True)
        
        #determine the new positive and negative set
        q_vals = uf.TDC_flex_c(SVM_train_labels == -1, SVM_train_labels == 1, c = 3/4, lam = 3/4)
        
        positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= alpha)
        
        negative_set_indxs = (SVM_train_labels == -1)
        
        SVM_train_features_iter = SVM_train_features.loc[positive_set_indxs | negative_set_indxs, :].copy()
        SVM_train_labels_iter = SVM_train_labels.loc[positive_set_indxs | negative_set_indxs].copy()
        
        #get actual power if we were to stop here
        real_labels = real_df['Label'].copy()
        real_df_test = real_df.drop(['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins'],axis = 1).copy()
        sds = real_df_test.apply(np.std, axis = 0)
        real_df_test = real_df_test.loc[:, abs(sds) > 1e-10]
        
        #scale non-binary features
        scale = StandardScaler()
        real_df_test.loc[:,~(real_df_test.columns.isin(["Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzN", "enzC", "rank", "n_o"]))] = scale.fit_transform(real_df_test.loc[:,~(real_df_test.columns.isin(["Charge1", "Charge2", "Charge3", "Charge4", "Charge5", "enzN", "enzC", "rank", "n_o"]))])
        
        #the new direction
        new_scores = grid.decision_function(real_df_test)
        new_idx = pd.Series(new_scores).sort_values(ascending = False).index
        new_labels = real_labels.loc[new_idx].reset_index(drop = True)
        
        q_val = uf.TDC_flex_c(new_labels == -1, new_labels == 1, c = 2/3, lam = 2/3)
        power_final = sum((q_val <= alpha) & (new_labels == 1))
        true_power[iterate] = power_final
        
        logging.info("Observed power: %s." %(power_final))
        sys.stderr.write("Observed power: %s.\n" %(power_final))
        
        logging.info("Trained power: %s." %(best_train_power))
        sys.stderr.write("Trained power: %s.\n" %(best_train_power))
        
        logging.info("Std trained power: %s." %(best_train_std))
        sys.stderr.write("Std trained power: %s.\n" %(best_train_std))
    
        
    #using the last new_idx to report the discoveries
    real_df = real_df.loc[new_idx].reset_index(drop = True)
    real_df_discoveries = real_df[(q_val <= alpha) & (new_labels == 1)]
        
    return(train_power, train_std, true_power, real_df_discoveries)

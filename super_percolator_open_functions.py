#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 19:13:19 2023

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
from sklearn.model_selection import RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
################################################################################


def get_bin_freq_pi0(df_all, precursor_bin_width=1.0005079/4):

    #get pi_0, bins, freq
    delta_mass_max = max(abs(df_all.ExpMass - df_all.CalcMass))
    breaks_p = np.arange(0, delta_mass_max + 2*precursor_bin_width,
                         precursor_bin_width) - precursor_bin_width/2
    breaks_n = list(reversed(-breaks_p))
    breaks = pd.Series(
        breaks_n[0:(len(breaks_n) - 1)] + list(breaks_p), name='bins')
    digitized = np.digitize(df_all.ExpMass - df_all.CalcMass, breaks)
    df_all['bins'] = digitized  # binning

    bin_freq = df_all['bins'].value_counts()  # getting bin frequencies
    bin_freq = bin_freq.reset_index()
    bin_freq.columns = ['bins', 'freq']
    df_all = df_all.merge(bin_freq, how='left', on='bins')

    #get pi_0
    pi_0s = df_all.groupby('bins').apply(lambda x: sum(x.Label == 1)/x.freq)
    pi_0s = pi_0s.reset_index()
    pi_0s.drop('level_1', axis=1, inplace=True)
    pi_0s.columns = ['bins', 'pi_0']
    pi_0s.drop_duplicates(inplace=True)
    df_all = df_all.merge(pi_0s, how='left', on='bins')

    return(df_all)
#########################################################################################################


def PSM_level(target_file, decoy_file, top=1):
    
    if type(decoy_file) == type(None):
        df = target_file
    elif type(decoy_file) == list:
        for i in range(len(decoy_file)):
            decoy_file[i]['filename'] = i
        decoy_file_combined = pd.concat(decoy_file)
        target_file['filename'] = 0
        df = pd.concat([target_file, decoy_file_combined])
    else:
        df = pd.concat([target_file, decoy_file])  # combine

    df = df.sample(frac=1).reset_index(drop=True)  # randomly shuffle
    
    df = df.sort_values(by='TailorScore', ascending=False).reset_index(
        drop=True)  # sort by score
    
    df['rank'] = df.groupby(["ScanNr", "ExpMass"])["TailorScore"].rank("first", ascending=False)

    df = df[df['rank'] <= top]  # get top PSMs for each scan
    
    df['rank'] = df['rank'].astype(int)

    df['SpecId'] = df.apply(lambda x: '_'.join(x.SpecId.split(
        '_')[:-1]) + '_' + str(x['rank']), axis=1)  # split specID

    return(df.reset_index(drop=True))
#########################################################################################################


def pseudo_PSM_level(df_all, df_extra_decoy, top, precursor_bin_width):

    df_all = df_all.loc[df_all['rank'] <= top, :].reset_index(drop=True).copy()

    df_all_temp = df_all.copy()  # copy to be reported separately

    # column values to take subset of df_extra_decoy
    keys = ['filename', 'ScanNr', 'rank']

    indx_all = df_all.set_index(keys).index  # index of df_all

    indx_extra_decoy = df_extra_decoy.set_index(
        keys).index  # index of df_extra_decoy

    df_extra_decoy = df_extra_decoy[indx_extra_decoy.isin(
        indx_all)]  # get subset of df_extra_decoy

    df_all['Label'] = 1  # create a pseudo target label

    # create combined dataframe
    df_combined = pd.concat([df_all, df_extra_decoy])

    df_combined = get_bin_freq_pi0(df_combined.copy(), precursor_bin_width)  # get extra features

    df_all_temp[['bins', 'freq', 'pi_0']
                ] = df_combined.loc[df_combined.Label == 1, ['bins', 'freq', 'pi_0']]

    df_combined = df_combined.sample(frac=1)  # randomly break ties

    df_combined = df_combined.sort_values(by='TailorScore', ascending=False).reset_index(
        drop=True)  # sort by score

    df_combined = df_combined.drop_duplicates(
        subset=['ScanNr', 'filename', 'rank', 'ExpMass'])  # doing pseudo psm level competition

    df_combined.reset_index(drop=True, inplace=True)

    return(df_combined, df_all_temp)
#########################################################################################################


def peptide_level(df_all, peptide_list_df, precursor_bin_width=1.0005079/4, keep_original = False, before = False, original_df = None):
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
    if 'n_o' in df_all.columns:
        df_all['rank'] = df_all['rank'].astype(int)
        df_all = df_all[((df_all['rank'] <= 2) & (df_all['n_o'] == 0)) | (
            (df_all['rank'] == 1) & (df_all['n_o'] == 1))]
        df_all.loc[df_all['rank'] == 1, 'rank'] = 0
        df_all.loc[df_all['rank'] == 2, 'rank'] = 1
    else:
        #taking top 1 PSMs
        sys.stderr.write("Taking top 1 PSMs. \n")
        logging.info("Taking top 1 PSMs.")
        df_all['rank'] = df_all['SpecId'].apply(
            lambda x: int(x[-1]))
        df_all = df_all[df_all['rank'] == 1].reset_index(drop = True)

    df_all = df_all.sample(frac=1)  # break ties randomly
    
    if 'SVM_score' in df_all.columns:
        df_all = df_all.sort_values(by='SVM_score', ascending=False).reset_index(
            drop=True)  # sort by score
    elif 'TailorScore' in df_all.columns:
        df_all = df_all.sort_values(by='TailorScore', ascending=False).reset_index(
            drop=True)  # sort by score
    else:
        df_all = df_all.sort_values(by='XCorr', ascending=False).reset_index(
            drop=True)  # sort by score

    # getting best score for each Peptide
    df_all = df_all.drop_duplicates(subset='Peptide')
    
    if type(peptide_list_df) == list:
        df_all['original_target'] = df_all['Peptide']
        for i in range(len(peptide_list_df)):
            df_all_sub = df_all[(df_all.Label == -1) & (df_all.filename == i)].copy()
            peptide_list_df[i].rename(
                columns={'target': 'original_target', 'decoy': 'Peptide'}, inplace=True)
            df_all_sub = df_all_sub.merge(
                peptide_list_df[i][['original_target', 'Peptide']], how='left', on='Peptide')
            df_all.loc[(df_all.Label == -1) & (df_all.filename == i),
                       'original_target'] = df_all_sub['original_target_y'].tolist()
    else:
        df_all_sub = df_all[df_all.Label == -1].copy()
        peptide_list_df.rename(
            columns={'target': 'original_target', 'decoy': 'Peptide'}, inplace=True)
        df_all_sub = df_all_sub.merge(
            peptide_list_df[['original_target', 'Peptide']], how='left', on='Peptide')
        df_all['original_target'] = df_all['Peptide']
        df_all.loc[df_all.Label == -1,
                   'original_target'] = df_all_sub['original_target'].tolist()
    
    df_all['original_target'] = df_all['original_target'].str.replace(
        "\\[|\\]|\\.|\\d+", "", regex=True)

    #adding both target-decoy scores so that they can be trained on as well
    if 'TailorScore' in df_all.columns:
        df_all = df_all.assign(min_tailor_score=df_all.groupby(
            'original_target').TailorScore.transform(lambda x: min(x) if min(x) != max(x) else 0))
    if 'XCorr' in df_all.columns:
        df_all = df_all.assign(min_xcorr_score=df_all.groupby(
            'original_target').XCorr.transform(lambda x: min(x) if min(x) != max(x) else 0))
    
    if before:
        df_all = pd.concat([df_all.loc[df_all.Label == 1], df_all.loc[df_all.Label == -1].drop_duplicates('original_target')])
        
        delta_mass_max = max(abs(df_all.ExpMass - df_all.CalcMass))
        breaks_p = np.arange(0, delta_mass_max + 2*precursor_bin_width,
                             precursor_bin_width) - precursor_bin_width/2
        breaks_n = list(reversed(-breaks_p))
        breaks = pd.Series(
            breaks_n[0:(len(breaks_n) - 1)] + list(breaks_p), name='bins')
        digitized = np.digitize(df_all.ExpMass - df_all.CalcMass, breaks)
        df_all['bins'] = digitized

        bin_freq = df_all['bins'].value_counts()  # getting bin frequencies
        bin_freq = bin_freq.reset_index()
        bin_freq.columns = ['bins', 'freq']
        df_all = df_all.merge(bin_freq, how='left', on='bins')
        
        if type(original_df) != type(None):
            original_df = original_df.merge(df_all[['SpecId', 'filename', 'ExpMass', 'freq', 'bins']], how = 'left', on = ['SpecId', 'filename', 'ExpMass'])
        
        if 'SVM_score' in df_all.columns:
            df_all = df_all.sort_values(by='SVM_score', ascending=False).reset_index(
                drop=True)  # sort by score
        elif 'TailorScore' in df_all.columns:
            df_all = df_all.sort_values(by='TailorScore', ascending=False).reset_index(
                drop=True)  # sort by score
        else:
            df_all = df_all.sort_values(by='XCorr', ascending=False).reset_index(
                drop=True)  # sort by score
        
        df_all = df_all.drop_duplicates(subset='original_target')
        if not keep_original:
            df_all.drop(['original_target', 'enzInt'], axis=1, inplace=True)
            if type(original_df) != type(None):
                original_df.drop(['original_target', 'enzInt'], axis=1, inplace=True, errors = 'ignore')
            
    else:
        
        df_all = df_all.drop_duplicates(subset='original_target')
        if not keep_original:
            df_all.drop(['original_target', 'enzInt'], axis=1, inplace=True, errors = 'ignore')
            if type(original_df) != type(None):
                original_df.drop(['original_target', 'enzInt'], axis=1, inplace=True, errors = 'ignore')
        # binning delta masses
        delta_mass_max = max(abs(df_all.ExpMass - df_all.CalcMass))
        breaks_p = np.arange(0, delta_mass_max + 2*precursor_bin_width,
                             precursor_bin_width) - precursor_bin_width/2
        breaks_n = list(reversed(-breaks_p))
        breaks = pd.Series(
            breaks_n[0:(len(breaks_n) - 1)] + list(breaks_p), name='bins')
        digitized = np.digitize(df_all.ExpMass - df_all.CalcMass, breaks)
        df_all['bins'] = digitized
    
        bin_freq = df_all['bins'].value_counts()  # getting bin frequencies
        bin_freq = bin_freq.reset_index()
        bin_freq.columns = ['bins', 'freq']
        df_all = df_all.merge(bin_freq, how='left', on='bins')
    
    if type(original_df) == type(None):
        return(df_all)
    else:
        return(df_all, original_df)
#########################################################################################################


def custom_accuracy(y, y_pred, alpha=0.01, p=0.5):
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
    y = [x for _, x in sorted(zip(y_pred, y), reverse=True)
         ]  # sorting y according to the predicted scores
    decoy_wins = [i == -1 for i in y]
    target_wins = [i == 1 for i in y]
    qvals = uf.TDC_flex_c(decoy_wins, target_wins, BC1=1,
                          c=1 - p/2, lam=1 - p/2)  # doing TDC
    disc = sum((qvals[i] <= alpha) and (y[i] == 1)
               for i in range(len(y)))  # getting discoveries
    return(disc)
#########################################################################################################


def train_cv(labels, df, folds=3, Cs=[0.1, 1, 10], kernel='linear', degree=2, alpha=0.01, p = 0.5):
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

    class_weight = [{-1: C_neg, 1: C_pos}
                    for C_neg in Cs for C_pos in Cs if C_neg >= C_pos]

    if kernel == 'linear':
        param_grid = dict(class_weight=class_weight)
    if kernel == 'rbf':
        gamma_range = np.logspace(-3, 5, 10)
        param_grid = dict(gamma=gamma_range, class_weight=class_weight)
    if kernel == 'poly':
        gamma_range = [0.1, 1, 10]
        coef0_range = [0.1, 1, 10]
        degree = [degree]
        param_grid = dict(gamma=gamma_range, degree=degree, coef0=coef0_range,
                          class_weight=class_weight)

    my_scorer = make_scorer(custom_accuracy, alpha=alpha, p = p,
                            greater_is_better=True, needs_threshold=True)

    if kernel == 'linear':
        grid = GridSearchCV(svm.SVC(), param_grid=param_grid,
                            cv=folds, scoring=my_scorer)
    else:
        grid = GridSearchCV(svm.SVC(), param_grid=param_grid,
                                  cv=folds, scoring=my_scorer)
    grid.fit(df, labels)
    while max(grid.cv_results_['mean_test_score']) == 0:
        alpha = alpha + 0.005
        my_scorer = make_scorer(custom_accuracy, alpha=alpha,
                                greater_is_better=True, needs_threshold=True)

        if kernel == 'linear':
            grid = GridSearchCV(svm.LinearSVC(max_iter = 1e5), param_grid=param_grid,
                                cv=folds, scoring=my_scorer)
        else:
            grid = RandomizedSearchCV(svm.SVC(), param_distributions=param_grid, n_iter=10,
                                      cv=folds, scoring=my_scorer)

        grid.fit(df, labels)

    return(grid)
#########################################################################################################


def train_lda_cv(labels, df):

    clf = LinearDiscriminantAnalysis()
    clf.fit(df, labels)

    return(clf)
#########################################################################################################

def do_scale(df_all, df_extra = None):
    #scale non-binary features
    scale = StandardScaler()
    if 'filename' in df_all.columns:
        df_all.loc[:, ~(df_all.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins', 'trained']))] = scale.fit_transform(
            df_all.loc[:, ~(df_all.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins', 'trained']))])
    else:
        df_all.loc[:, ~(df_all.columns.isin(['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins', 'trained']))] = scale.fit_transform(
            df_all.loc[:, ~(df_all.columns.isin(['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins', 'trained']))])
    
    if type(df_extra) == type(None):
        return(df_all)
    else:
        if 'filename' in df_extra.columns:
            df_extra.loc[:, ~(df_extra.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins', 'trained']))] = scale.transform(
                df_extra.loc[:, ~(df_extra.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins', 'trained']))])
        else:
            df_extra.loc[:, ~(df_extra.columns.isin(['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins', 'trained']))] = scale.transform(
                df_extra.loc[:, ~(df_extra.columns.isin(['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins', 'trained']))])
        return(df_all, df_extra)

################################################################################################

def do_svm(df_all, train_all, folds=3, Cs=[0.01, 0.1, 1, 10], total_iter=10, p=0.5, kernel='linear', alpha=0.01, train_alpha=0.01, degree=None, remove=None, top_positive=True, mult = 1):

    train_alpha_init = train_alpha
    
    train_targets = df_all[~(df_all.index.isin(train_all.index))].copy() #ensure the indices of train_all belong to df_all
    train_targets.loc[:, 'Label'] = 1

    train_df = pd.concat([train_all, train_targets]).reset_index(drop=True)
    train_df = train_df.sample(frac=1).reset_index(drop=True)
    
    if type(remove) == list:
        train_df.drop(remove, axis=1, inplace=True, errors = 'ignore')
    if 'SVM_score' in train_df.columns:
        train_df = train_df.sort_values(
            by = 'SVM_score', ascending=False).reset_index(drop = True)
    elif 'TailorScore' in train_df.columns:
        train_df = train_df.sort_values(
            by='TailorScore', ascending=False).reset_index(drop=True)
    else:
        train_df = train_df.sort_values(
            by='XCorr', ascending=False).reset_index(drop=True)

    #real df
    real_df = df_all[~(df_all.index.isin(train_all.index))
                     ].copy().reset_index(drop=True)
    if type(remove) == list:
        real_df.drop(remove, axis=1, inplace=True, errors = 'ignore')

    #Preprocess dataframe
    SVM_train_labels = train_df['Label'].copy()
    if 'filename' in train_df.columns:
        SVM_train_features = train_df.drop(
            ['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
    else:
        SVM_train_features = train_df.drop(
            ['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
        
    #Get rid of colinear features
    sds = SVM_train_features.apply(np.std, axis = 0)
    SVM_train_features = SVM_train_features[SVM_train_features.columns[sds != 0]]
    
    positive_set_indxs = [0]
    
    while sum(positive_set_indxs) == 0:
        train_alpha = train_alpha + 0.005
        #getting initial positive and negative set
        q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                               SVM_train_labels == 1, c=(mult*(1 - p) + 1)/(mult + 1), lam=(mult*(1 - p) + 1)/(mult + 1))
        if top_positive:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId)) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
        else:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId))
        negative_set_indxs = (SVM_train_labels == -1)
        
    train_alpha = train_alpha_init
    
    SVM_train_features_iter = SVM_train_features.loc[positive_set_indxs | negative_set_indxs, :].reset_index(
        drop=True).copy()
    SVM_train_labels_iter = SVM_train_labels.loc[positive_set_indxs | negative_set_indxs].reset_index(
        drop=True).copy()

    SVM_train_features_iter = SVM_train_features_iter.sample(frac=1)
    SVM_train_labels_iter = SVM_train_labels_iter.loc[SVM_train_features_iter.index]

    SVM_train_features_iter.reset_index(drop=True, inplace=True)
    SVM_train_labels_iter.reset_index(drop=True, inplace=True)

    train_power, train_std, true_power = [
        0]*total_iter, [0]*total_iter, [0]*total_iter

    logging.info("Conducting iterative SVM.")
    sys.stderr.write("Conducting iterative SVM.\n")

    for iterate in range(total_iter):
        logging.info("iteration: %s." % (iterate))
        sys.stderr.write("iteration: %s.\n" % (iterate))
        #determining best direction with cross validation for parameter selection
        grid = train_cv(SVM_train_labels_iter, SVM_train_features_iter,
                        folds=folds, Cs=Cs, kernel=kernel, degree=degree, alpha=train_alpha, p=p)
        best_train_power = max(grid.cv_results_['mean_test_score'])
        best_train_std = max(grid.cv_results_['std_test_score'])
        train_power[iterate] = best_train_power
        train_std[iterate] = best_train_std

        #the new direction
        new_scores = grid.decision_function(SVM_train_features)

        new_idx = pd.Series(new_scores).sort_values(ascending=False).index

        SVM_train_features = SVM_train_features.loc[new_idx].reset_index(
            drop=True)

        SVM_train_labels = SVM_train_labels.loc[new_idx].reset_index(drop=True)

        #determine the new positive and negative set
        q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                               SVM_train_labels == 1, c=(mult*(1 - p) + 1)/(mult + 1), lam=(mult*(1 - p) + 1)/(mult + 1))
        
        if top_positive:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId)) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
        else:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId))
        negative_set_indxs = (SVM_train_labels == -1)
        
        while sum(positive_set_indxs) == 0:
            train_alpha = train_alpha + 0.005
            #getting initial positive and negative set
            q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                                   SVM_train_labels == 1, c=(mult*(1 - p) + 1)/(mult + 1), lam=(mult*(1 - p) + 1)/(mult + 1))
            if top_positive:
                positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId)) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
            else:
                positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId))
            negative_set_indxs = (SVM_train_labels == -1)
        
        train_alpha = train_alpha_init
        
        SVM_train_features_iter = SVM_train_features.loc[positive_set_indxs | negative_set_indxs, :].reset_index(
            drop=True).copy()
        SVM_train_labels_iter = SVM_train_labels.loc[positive_set_indxs | negative_set_indxs].reset_index(
            drop=True).copy()

        SVM_train_features_iter = SVM_train_features_iter.sample(frac=1)
        SVM_train_labels_iter = SVM_train_labels_iter.loc[SVM_train_features_iter.index]

        SVM_train_features_iter.reset_index(drop=True, inplace=True)
        SVM_train_labels_iter.reset_index(drop=True, inplace=True)

        #get actual power if we were to stop here
        real_labels = real_df['Label'].copy()
        if 'filename' in real_df.columns:
            real_df_test = real_df.drop(
                ['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
        else:
            real_df_test = real_df.drop(
                ['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
            
        #Get rid of redundant features
        real_df_test = real_df_test[real_df_test.columns[real_df_test.columns.isin(SVM_train_features.columns)]]
    
        #the new direction
        new_scores = grid.decision_function(real_df_test)
        new_idx = pd.Series(new_scores).sort_values(ascending=False).index
        new_labels = real_labels.loc[new_idx].reset_index(drop=True)

        q_val = uf.TDC_flex_c(
            new_labels == -1, new_labels == 1, c=1/(mult*(1 - p) + 1), lam=1/(mult*(1 - p) + 1))
        power_final = sum((q_val <= alpha) & (new_labels == 1))
        true_power[iterate] = power_final

        logging.info("Observed power: %s." % (power_final))
        sys.stderr.write("Observed power: %s.\n" % (power_final))

        logging.info("Trained power: %s." % (best_train_power))
        sys.stderr.write("Trained power: %s.\n" % (best_train_power))

        logging.info("Std trained power: %s." % (best_train_std))
        sys.stderr.write("Std trained power: %s.\n" % (best_train_std))
        
        

    #using the last new_idx to report the discoveries
    real_df['SVM_score'] = new_scores
    
    if type(remove) == list:
        train_all_test = train_all.drop(remove, axis=1, errors = 'ignore')

    if 'filename' in train_all_test.columns:
        train_all_test = train_all_test.drop(
            ['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
    else:
        train_all_test = train_all_test.drop(
            ['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
        
    #Get rid of redundant features
    train_all_test = train_all_test[train_all_test.columns[train_all_test.columns.isin(SVM_train_features.columns)]]
    
    new_scores = grid.decision_function(train_all_test)
    
    train_all['SVM_score'] = new_scores

    return(train_power, train_std, true_power, real_df, train_all)


##########


def do_lda(df_all, train_decoys, total_iter=10, p=0.5, alpha=0.01, train_alpha=0.01, remove=None, top_positive=True, qda=False):

    train_targets = df_all[~(df_all.index.isin(train_decoys.index))].copy()
    train_targets.loc[:, 'Label'] = 1

    train_df = pd.concat([train_decoys, train_targets]).reset_index(drop=True)
    train_df = train_df.sample(frac=1).reset_index(drop=True)
    
    if type(remove) == list:
        train_df.drop(remove, axis=1, inplace=True, errors = 'ignore')
    if 'SVM_score' in train_df.columns:
        train_df = train_df.sort_values(
            by = 'SVM_score', ascending=False).reset_index(drop = True)
    elif 'TailorScore' in train_df.columns:
        train_df = train_df.sort_values(
            by='TailorScore', ascending=False).reset_index(drop=True)
    else:
        train_df = train_df.sort_values(
            by='XCorr', ascending=False).reset_index(drop=True)

    #real df
    real_df = df_all[~(df_all.index.isin(train_decoys.index))
                     ].copy().reset_index(drop=True)
    if type(remove) == list:
        real_df.drop(remove, axis=1, inplace=True, errors = 'ignore')

    #Preprocess dataframe
    SVM_train_labels = train_df['Label'].copy()
    if 'filename' in train_df.columns:
        SVM_train_features = train_df.drop(
            ['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
    else:
        SVM_train_features = train_df.drop(
            ['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
        
    #Get rid of colinear features
    sds = SVM_train_features.apply(np.std, axis = 0)
    SVM_train_features = SVM_train_features[SVM_train_features.columns[sds != 0]]

    #getting initial positive and negative set
    q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                           SVM_train_labels == 1, c=1 - p/2, lam=1 - p/2)
    if top_positive:
        positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId)) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
    else:
        positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId))
    negative_set_indxs = (SVM_train_labels == -1)
    
    while sum(positive_set_indxs) == 0:
        train_alpha = train_alpha + 0.005
        #getting initial positive and negative set
        q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                               SVM_train_labels == 1, c=1 - p/2, lam=1 - p/2)
        if top_positive:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId)) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
        else:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId))
        negative_set_indxs = (SVM_train_labels == -1)

    SVM_train_features_iter = SVM_train_features.loc[positive_set_indxs | negative_set_indxs, :].reset_index(
        drop=True).copy()
    SVM_train_labels_iter = SVM_train_labels.loc[positive_set_indxs | negative_set_indxs].reset_index(
        drop=True).copy()

    SVM_train_features_iter = SVM_train_features_iter.sample(frac=1)
    SVM_train_labels_iter = SVM_train_labels_iter.loc[SVM_train_features_iter.index]

    SVM_train_features_iter.reset_index(drop=True, inplace=True)
    SVM_train_labels_iter.reset_index(drop=True, inplace=True)

    train_power, train_std, true_power = [
        0]*total_iter, [0]*total_iter, [0]*total_iter

    logging.info("Conducting iterative lda.")
    sys.stderr.write("Conducting iterative lda.\n")

    for iterate in range(total_iter):
        logging.info("iteration: %s." % (iterate))
        sys.stderr.write("iteration: %s.\n" % (iterate))
        #determining best direction with cross validation for parameter selection
        grid = train_lda_cv(SVM_train_labels_iter, SVM_train_features_iter)
        
        #best_train_power = max(grid.cv_results_['mean_test_score'])
        #best_train_std = max(grid.cv_results_['std_test_score'])
        #train_power[iterate] = best_train_power
        #train_std[iterate] = best_train_std

        #the new direction
        new_scores = grid.decision_function(SVM_train_features)

        new_idx = pd.Series(new_scores).sort_values(ascending=False).index

        SVM_train_features = SVM_train_features.loc[new_idx].reset_index(
            drop=True)

        SVM_train_labels = SVM_train_labels.loc[new_idx].reset_index(drop=True)

        #determine the new positive and negative set
        q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                               SVM_train_labels == 1, c=1 - p/2, lam=1 - p/2)
        
        if top_positive:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= alpha) & (train_df.SpecId.isin(real_df.SpecId)) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
        else:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= alpha) & (train_df.SpecId.isin(real_df.SpecId))
        negative_set_indxs = (SVM_train_labels == -1)
        
        while sum(positive_set_indxs) == 0:
            train_alpha = train_alpha + 0.005
            #getting initial positive and negative set
            q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                                   SVM_train_labels == 1, c=1 - p/2, lam=1 - p/2)
            if top_positive:
                positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId)) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
            else:
                positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId))
            negative_set_indxs = (SVM_train_labels == -1)

        SVM_train_features_iter = SVM_train_features.loc[positive_set_indxs | negative_set_indxs, :].reset_index(
            drop=True).copy()
        SVM_train_labels_iter = SVM_train_labels.loc[positive_set_indxs | negative_set_indxs].reset_index(
            drop=True).copy()

        SVM_train_features_iter = SVM_train_features_iter.sample(frac=1)
        SVM_train_labels_iter = SVM_train_labels_iter.loc[SVM_train_features_iter.index]

        SVM_train_features_iter.reset_index(drop=True, inplace=True)
        SVM_train_labels_iter.reset_index(drop=True, inplace=True)

        #get actual power if we were to stop here
        real_labels = real_df['Label'].copy()
        if 'filename' in real_df.columns:
            real_df_test = real_df.drop(
                ['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
        else:
            real_df_test = real_df.drop(
                ['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
            
        #Get rid of colinear features
        sds = real_df_test.apply(np.std, axis = 0)
        real_df_test = real_df_test[real_df_test.columns[sds != 0]]

        #the new direction
        new_scores = grid.decision_function(real_df_test)
        new_idx = pd.Series(new_scores).sort_values(ascending=False).index
        new_labels = real_labels.loc[new_idx].reset_index(drop=True)

        q_val = uf.TDC_flex_c(
            new_labels == -1, new_labels == 1, c=1/(2 - p), lam=1/(2 - p))
        power_final = sum((q_val <= alpha) & (new_labels == 1))
        true_power[iterate] = power_final

        logging.info("Observed power: %s." % (power_final))
        sys.stderr.write("Observed power: %s.\n" % (power_final))

        #logging.info("Trained power: %s." % (best_train_power))
        #sys.stderr.write("Trained power: %s.\n" % (best_train_power))

        #logging.info("Std trained power: %s." % (best_train_std))
        #sys.stderr.write("Std trained power: %s.\n" % (best_train_std))

    #using the last new_idx to report the discoveries
    real_df['SVM_score'] = new_scores
    
    if type(remove) == list:
        train_decoys_test = train_decoys.drop(remove, axis=1, errors = 'ignore')

    if 'filename' in train_decoys_test.columns:
        train_decoys_test = train_decoys_test.drop(
            ['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
    else:
        train_decoys_test = train_decoys_test.drop(
            ['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
        
    #Get rid of colinear features
    sds = train_decoys_test.apply(np.std, axis = 0)
    train_decoys_test = train_decoys_test[train_decoys_test.columns[sds != 0]]
    
    new_scores = grid.decision_function(train_decoys_test)
    
    train_decoys['SVM_score'] = new_scores

    return(real_df, train_decoys)
    

################################################################################################

def do_svm_extra(df_pseudo, df_all, folds=3, Cs=[0.01, 0.1, 1, 10], total_iter=10, kernel='linear', alpha=0.01, train_alpha=0.01, degree=None, remove=None, top_positive=True):
    train_alpha_init = train_alpha
    
    train_df = df_pseudo.copy().sample(frac=1).reset_index(drop=True)
    
    if type(remove) == list:
        train_df.drop(remove, axis=1, inplace=True, errors = 'ignore')
    if 'SVM_score' in train_df.columns:
        train_df = train_df.sort_values(
            by = 'SVM_score', ascending=False).reset_index(drop = True)
    elif 'TailorScore' in train_df.columns:
        train_df = train_df.sort_values(
            by='TailorScore', ascending=False).reset_index(drop=True)
    else:
        train_df = train_df.sort_values(
            by='XCorr', ascending=False).reset_index(drop=True)

    #real df
    real_df = df_all.copy().reset_index(drop=True)
    if type(remove) == list:
        real_df.drop(remove, axis=1, inplace=True, errors = 'ignore')

    #Preprocess dataframe
    SVM_train_labels = train_df['Label'].copy()
    if 'filename' in train_df.columns:
        SVM_train_features = train_df.drop(
            ['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
    else:
        SVM_train_features = train_df.drop(
            ['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
        
    #Get rid of colinear features
    sds = SVM_train_features.apply(np.std, axis = 0)
    SVM_train_features = SVM_train_features[SVM_train_features.columns[sds != 0]]
    
    positive_set_indxs = [0]
    
    while sum(positive_set_indxs) == 0:
        
        #getting initial positive and negative set
        q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                               SVM_train_labels == 1, c=2/3, lam=2/3)
        if top_positive:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
        else:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha)
        negative_set_indxs = (SVM_train_labels == -1)
        
        train_alpha = train_alpha + 0.005
        
    train_alpha = train_alpha_init
    
    SVM_train_features_iter = SVM_train_features.loc[positive_set_indxs | negative_set_indxs, :].reset_index(
        drop=True).copy()
    SVM_train_labels_iter = SVM_train_labels.loc[positive_set_indxs | negative_set_indxs].reset_index(
        drop=True).copy()

    SVM_train_features_iter = SVM_train_features_iter.sample(frac=1)
    SVM_train_labels_iter = SVM_train_labels_iter.loc[SVM_train_features_iter.index]

    SVM_train_features_iter.reset_index(drop=True, inplace=True)
    SVM_train_labels_iter.reset_index(drop=True, inplace=True)

    train_power, train_std, true_power = [
        0]*total_iter, [0]*total_iter, [0]*total_iter

    logging.info("Conducting iterative SVM.")
    sys.stderr.write("Conducting iterative SVM.\n")

    for iterate in range(total_iter):
        logging.info("iteration: %s." % (iterate))
        sys.stderr.write("iteration: %s.\n" % (iterate))
        #determining best direction with cross validation for parameter selection
        grid = train_cv(SVM_train_labels_iter, SVM_train_features_iter,
                        folds=folds, Cs=Cs, kernel=kernel, degree=degree, alpha=train_alpha, p=2/3)
        best_train_power = max(grid.cv_results_['mean_test_score'])
        best_train_std = max(grid.cv_results_['std_test_score'])
        train_power[iterate] = best_train_power
        train_std[iterate] = best_train_std

        #the new direction
        new_scores = grid.decision_function(SVM_train_features)

        new_idx = pd.Series(new_scores).sort_values(ascending=False).index

        SVM_train_features = SVM_train_features.loc[new_idx].reset_index(
            drop=True)

        SVM_train_labels = SVM_train_labels.loc[new_idx].reset_index(drop=True)

        #determine the new positive and negative set
        q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                               SVM_train_labels == 1, c=2/3, lam=2/3)
        
        if top_positive:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
        else:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha)
        negative_set_indxs = (SVM_train_labels == -1)
        
        while sum(positive_set_indxs) == 0:
            train_alpha = train_alpha + 0.005
            #getting initial positive and negative set
            q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                                   SVM_train_labels == 1, c=2/3, lam=2/3)
            if top_positive:
                positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
            else:
                positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha)
            negative_set_indxs = (SVM_train_labels == -1)
        
        train_alpha = train_alpha_init
        
        SVM_train_features_iter = SVM_train_features.loc[positive_set_indxs | negative_set_indxs, :].reset_index(
            drop=True).copy()
        SVM_train_labels_iter = SVM_train_labels.loc[positive_set_indxs | negative_set_indxs].reset_index(
            drop=True).copy()

        SVM_train_features_iter = SVM_train_features_iter.sample(frac=1)
        SVM_train_labels_iter = SVM_train_labels_iter.loc[SVM_train_features_iter.index]

        SVM_train_features_iter.reset_index(drop=True, inplace=True)
        SVM_train_labels_iter.reset_index(drop=True, inplace=True)

        #get actual power if we were to stop here
        real_labels = real_df['Label'].copy()
        if 'filename' in real_df.columns:
            real_df_test = real_df.drop(
                ['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
        else:
            real_df_test = real_df.drop(
                ['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
            
        #Get rid of redundant features
        real_df_test = real_df_test[real_df_test.columns[real_df_test.columns.isin(SVM_train_features.columns)]]
    
        #the new direction
        new_scores = grid.decision_function(real_df_test)
        new_idx = pd.Series(new_scores).sort_values(ascending=False).index
        new_labels = real_labels.loc[new_idx].reset_index(drop=True)

        q_val = uf.TDC_flex_c(
            new_labels == -1, new_labels == 1, c=1/2, lam=1/2)
        power_final = sum((q_val <= alpha) & (new_labels == 1))
        true_power[iterate] = power_final

        logging.info("Observed power: %s." % (power_final))
        sys.stderr.write("Observed power: %s.\n" % (power_final))

        logging.info("Trained power: %s." % (best_train_power))
        sys.stderr.write("Trained power: %s.\n" % (best_train_power))

        logging.info("Std trained power: %s." % (best_train_std))
        sys.stderr.write("Std trained power: %s.\n" % (best_train_std))
        
        

    #using the last new_idx to report the discoveries
    real_df['SVM_score'] = new_scores

    return(train_power, train_std, true_power, real_df)

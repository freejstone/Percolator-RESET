#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 12:47:06 2023

@author: jfre0619
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


def get_rank(df, top):

    df = df.sample(frac=1)

    df['rank'] = df.groupby(["filename", "ScanNr", "ExpMass"])["XCorr"].rank(
        "first", ascending=False).astype(int)  # get ranking

    df = df[df['rank'] <= top]  # get top PSMs for each scan

    df.reset_index(drop=True, inplace=True)

    return(df)

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

################################################################################


def pseudo_PSM_level(df_all, df_extra_decoy, top, precursor_bin_width):

    df_all = df_all.loc[df_all['rank'] <= top, :].reset_index(drop=True).copy()

    df_all_temp = df_all.copy()  # copy to be reported separately

    # column values to take subset of df_extra_decoy
    keys = ['filename', 'ScanNr', 'n_o', 'rank']

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
        subset=['ScanNr', 'filename', 'n_o', 'rank', 'ExpMass'])  # doing pseudo psm level competition

    df_combined.reset_index(drop=True, inplace=True)

    return(df_combined, df_all_temp)


################################################################################
def peptide_level(df_all, peptide_list_df):
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

    df_all = df_all.sample(frac=1)  # break ties randomly

    df_all = df_all.sort_values(by='SVM_scores', ascending=False).reset_index(
        drop=True)  # sort by score

    df_all_sub = df_all[df_all.Label == -1].copy()

    if 'original_target' in df_all_sub.columns:
        df_all_sub.drop('original_target', axis=1, inplace=True)

    peptide_list_df.rename(
        columns={'target': 'original_target', 'decoy': 'Peptide'}, inplace=True)
    df_all_sub = df_all_sub.merge(
        peptide_list_df[['original_target', 'Peptide']], how='left', on='Peptide')

    if 'original_target' not in df_all.columns:
        df_all['original_target'] = df_all['Peptide']
    df_all.loc[df_all.Label == -1,
               'original_target'] = df_all_sub['original_target'].tolist()

    df_all['original_target'] = df_all['original_target'].str.replace(
        "\\[|\\]|\\.|\\d+", "", regex=True)

    df_all = df_all.drop_duplicates(subset='original_target')

    return(df_all)

################################################################################


def PSM_level(target_file, decoy_file, top=1):
    df = pd.concat([target_file, decoy_file])  # combine

    df = df.sample(frac=1).reset_index(drop=True)  # randomly shuffle

    df['rank'] = df.groupby(["filename", "ScanNr", "ExpMass"])["XCorr"].rank(
        "first", ascending=False).astype(int)  # get ranking

    df = df[df['rank'] <= top]  # get top PSMs for each scan

    df['SpecId'] = df.apply(lambda x: '_'.join(x.SpecId.split(
        '_')[:-1]) + '_' + str(x['rank']), axis=1)  # split specID

    return(df.reset_index(drop=True))

################################################################################


def custom_accuracy(y, y_pred, alpha=0.01):
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
                          c=3/4, lam=3/4)  # doing TDC
    disc = sum((qvals[i] <= alpha) and (y[i] == 1)
               for i in range(len(y)))  # getting discoveries
    return(disc)
#########################################################################################################


def train_cv(labels, df, folds=3, Cs=[0.1, 1, 10], kernel='linear', degree=2, alpha=0.01):
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

    my_scorer = make_scorer(custom_accuracy, alpha=alpha,
                            greater_is_better=True, needs_threshold=True)

    if kernel == 'linear':
        grid = GridSearchCV(svm.SVC(), param_grid=param_grid,
                            cv=folds, scoring=my_scorer)
    else:
        grid = RandomizedSearchCV(svm.SVC(), param_distributions=param_grid, n_iter=10,
                                  cv=folds, scoring=my_scorer)
    grid.fit(df, labels)
    while max(grid.cv_results_['mean_test_score']) == 0:
        alpha = alpha + 0.005
        my_scorer = make_scorer(custom_accuracy, alpha=alpha,
                                greater_is_better=True, needs_threshold=True)

        if kernel == 'linear':
            grid = GridSearchCV(svm.SVC(), param_grid=param_grid,
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


def do_iterative_svm_cv(df_all, real_df, folds=3, Cs=[0.1, 1, 10], total_iter=10, kernel='linear', alpha=0.01, train_alpha=0.01, degree=None, remove=None, top_positive=True):

    train_df = df_all.copy()
    train_df = train_df.sample(frac=1)
    if type(remove) == list:
        train_df.drop(remove, axis=1, inplace=True)
    if 'TailorScore' in train_df.columns:
        train_df = train_df.sort_values(
            by='TailorScore', ascending=False).reset_index(drop=True)
    else:
        train_df = train_df.sort_values(
            by='XCorr', ascending=False).reset_index(drop=True)

    #scale non-character features
    scale = StandardScaler()
    if 'filename' in train_df.columns:
        train_df.loc[:, ~(train_df.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins']))] = scale.fit_transform(
            train_df.loc[:, ~(train_df.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins']))])
    else:
        train_df.loc[:, ~(train_df.columns.isin(['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins']))] = scale.fit_transform(
            train_df.loc[:, ~(train_df.columns.isin(['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins']))])

    #real df
    if type(remove) == list:
        real_df.drop(remove, axis=1, inplace=True)

    #scale real_df according to the fitted scale
    real_df.loc[:, ~(real_df.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins']))] = scale.transform(
        real_df.loc[:, ~(real_df.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins']))])

    #Preprocess dataframe
    SVM_train_labels = train_df['Label'].copy()
    if 'filename' in train_df.columns:
        SVM_train_features = train_df.drop(
            ['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()
    else:
        SVM_train_features = train_df.drop(
            ['SpecId', 'Label', 'ScanNr', 'Peptide', 'Proteins'], axis=1).copy()

    #Get rid of colinear features
    sds = SVM_train_features.apply(np.std, axis=0)
    SVM_train_features = SVM_train_features[SVM_train_features.columns[sds != 0]]

    #getting initial positive and negative set
    q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                           SVM_train_labels == 1, c=2/3, lam=2/3)
    if top_positive:
        positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (
            SVM_train_features['rank'] == min(SVM_train_features['rank']))
    else:
        positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha)
    negative_set_indxs = (SVM_train_labels == -1)

    while sum(positive_set_indxs) == 0:
        train_alpha = train_alpha + 0.005
        #getting initial positive and negative set
        q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                               SVM_train_labels == 1, c=2/3, lam=2/3)
        if top_positive:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (
                SVM_train_features['rank'] == min(SVM_train_features['rank']))
        else:
            positive_set_indxs = (SVM_train_labels == 1) & (
                q_vals <= train_alpha)
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

    logging.info("Conducting iterative SVM.")
    sys.stderr.write("Conducting iterative SVM.\n")

    for iterate in range(total_iter):
        logging.info("iteration: %s." % (iterate))
        sys.stderr.write("iteration: %s.\n" % (iterate))
        #determining best direction with cross validation for parameter selection
        grid = train_cv(SVM_train_labels_iter, SVM_train_features_iter,
                        folds=folds, Cs=Cs, kernel=kernel, degree=degree, alpha=train_alpha)
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
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= alpha) & (
                SVM_train_features['rank'] == min(SVM_train_features['rank']))
        else:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= alpha)
        negative_set_indxs = (SVM_train_labels == -1)

        while sum(positive_set_indxs) == 0:
            train_alpha = train_alpha + 0.005
            #getting initial positive and negative set
            q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                                   SVM_train_labels == 1, c=2/3, lam=2/3)
            if top_positive:
                positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (
                    SVM_train_features['rank'] == min(SVM_train_features['rank']))
            else:
                positive_set_indxs = (SVM_train_labels == 1) & (
                    q_vals <= train_alpha)
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
        sds = real_df_test.apply(np.std, axis=0)
        real_df_test = real_df_test[real_df_test.columns[sds != 0]]

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
    real_df['SVM_scores'] = new_scores

    return(train_power, train_std, true_power, real_df)
#########################################################################################################

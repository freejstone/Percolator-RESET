#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 19:13:19 2023

@author: jackfreestone

super_percolator main functions
"""

import numpy as np
import pandas as pd
import utility_functions as uf
import logging
import sys
from sklearn import svm
from sklearn.metrics import make_scorer
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.preprocessing import StandardScaler
################################################################################


def PSM_level(target_file, decoy_file, top=1):
    '''
    

    Parameters
    ----------
    target_file : Pandas Dataframe
        Target search file when concat=F.
    decoy_file : Pandas Dataframe
        Decoy search file when concat=F.
    top : Int, optional
        DESCRIPTION. The default is 1. The number of top PSMs to consider from
        either search file.

    Returns
    -------
    The top 1 PSMs searched against the concatenated databases. Also adjusts
    some features so that they are measured with respect to the concantenated
    databases. Removes enzInt feature.

    '''
    if decoy_file is None:
        df = target_file
    elif isinstance(decoy_file, list):
        for i in range(len(decoy_file)):
            decoy_file[i]['filename'] = i
        decoy_file_combined = pd.concat(decoy_file)
        target_file['filename'] = 0
        df = pd.concat([target_file, decoy_file_combined])
    else:
        df = pd.concat([decoy_file, target_file])  # combine

    df = df.sample(frac=1).reset_index(drop=True)  # randomly shuffle

    df['Charge'] = df['SpecId'].apply(
        lambda x: int(x[-3]))

    #Fix up lnNumSP
    logging.info(
        "Fixing up lnNumSP so that it is with respect to the target database.")
    sys.stderr.write(
        "Fixing up lnNumSP so that it is with respect to the target database.\n")
    values_x = df.loc[df['Label'] == 1].groupby(
        ["ScanNr", "ExpMass", "Charge"])['lnNumSP'].first()
   
    df['lnNumSP'] = df.apply(lambda x: values_x[(
        x['ScanNr'], x['ExpMass'], x['Charge'])], axis=1)

    df = df.sort_values(by='XCorr', ascending=False).reset_index(
        drop=True)  # sort by score

    df.drop_duplicates(["ScanNr", "ExpMass", "Charge"])

    df['rank'] = df.groupby(["ScanNr", "ExpMass", "Charge"])[
        "XCorr"].rank("first", ascending=False)

    df['rank'] = df['rank'].astype(int)

    df = df[df['rank'] <= 5]  # get top PSMs for each scan

    df['SpecId'] = df.apply(lambda x: '_'.join(x.SpecId.split(
        '_')[:-1]) + '_' + str(x['rank']), axis=1)  # split specID

    #Fix up deltlCn
    logging.info(
        "Fixing up deltlCn so that it is with respect to the concatenated database.")
    sys.stderr.write(
        "Fixing up deltlCn so that it is with respect to the concatenated database.\n")
    LastValue = df.groupby(["ScanNr", "ExpMass", "Charge"])[
        'XCorr'].transform('last')
    df['deltLCn'] = (df['XCorr'] - LastValue)/np.maximum(df['XCorr'], 1)

    #Fix up deltCn
    logging.info(
        "Fixing up deltCn so that it is with respect to the concatenated database.")
    sys.stderr.write(
        "Fixing up deltCn so that it is with respect to the concatenated database.\n")
    df['SubsequentValue'] = df.groupby(["ScanNr", "ExpMass", "Charge"])[
        'XCorr'].shift(-1)
    df['deltCn'] = df['XCorr'].sub(
        df['SubsequentValue'])/np.maximum(df['XCorr'], 1)
    
    df.loc[df['deltCn'].isna(), 'deltCn'] = 0

    df.drop(['Charge', 'SubsequentValue'], axis=1, inplace=True)

    df = df[df['rank'] <= top]
    
    df.drop(['rank'], axis=1, inplace=True)

    return(df.reset_index(drop=True))
#########################################################################################################


def peptide_level(df_all, peptide_list_df, pair, initial_dir):
    '''

    Parameters
    ----------
    df_all : Pandas Dataframe
        Concatenated search file.
    peptide_list_df : Pandas Dataframe
        Tide-index target-decoy peptide pairs.
    Returns
    -------
    PSMs that remain after target-decoy peptide level competition.

    '''

    #select only top 1 PSMs
    sys.stderr.write("Taking top 1 PSMs. \n")
    logging.info("Taking top 1 PSMs.")
    df_all['rank'] = df_all['SpecId'].apply(
        lambda x: int(x[-1]))
    df_all = df_all[df_all['rank'] == 1].reset_index(drop=True)

    df_all = df_all.sample(frac=1).reset_index(
        drop=True)  # break ties randomly
    
    get_columns = df_all.columns[df_all.columns.str.contains(initial_dir, case=False)]
    
    if len(get_columns) > 0:
        initial_dir_case_insensitive = df_all.columns[df_all.columns.str.contains(initial_dir, case=False)][0]
    else:
        sys.exit("--initial_dir %s not detected. \n" %(initial_dir))
        
    
    df_all = df_all.sort_values(by=initial_dir_case_insensitive, ascending=False).reset_index(
        drop=True)  # sort by score
        

    # getting best score for each Peptide
    df_all = df_all.drop_duplicates(subset='Peptide')

    sys.stderr.write("Doing peptide-stem level competition. \n")
    logging.info("Doing peptide-stem level competition.")
    
    if not pair:
        df_all['original_target'] = df_all['Peptide']
    
    else:
        if isinstance(peptide_list_df, list):
            df_all['original_target'] = df_all['Peptide']
            for i in range(len(peptide_list_df)):
                df_all_sub = df_all[(df_all.Label == -1) &
                                    (df_all.filename == i)].copy()
                peptide_list_df[i].rename(
                    columns={'target': 'original_target', 'decoy': 'Peptide'}, inplace=True)
                df_all_sub = df_all_sub.merge(
                    peptide_list_df[i][['original_target', 'Peptide']], how='left', on='Peptide')
                
                if any(df_all_sub.original_target_y.isna()):
                    sys.exit("Some peptides in the search file do not have a pair in the peptide list. E.g. %s \n" %(df_all_sub.Peptide[df_all_sub.original_target_y.isna()].values[0])) 
                    
                df_all.loc[(df_all.Label == -1) & (df_all.filename == i),
                           'original_target'] = df_all_sub['original_target_y'].tolist()
        else:
            df_all_sub = df_all[df_all.Label == -1].copy()
            peptide_list_df.rename(
                columns={'target': 'original_target', 'decoy': 'Peptide'}, inplace=True)
            df_all_sub = df_all_sub.merge(
                peptide_list_df[['original_target', 'Peptide']], how='left', on='Peptide')
            
            if any(df_all_sub.original_target.isna()):
                sys.exit("Some peptides in the search file do not have a pair in the peptide list. E.g. %s \n" %(df_all_sub.Peptide[df_all_sub.original_target.isna()].values[0])) 
                
            df_all['original_target'] = df_all['Peptide']
            df_all.loc[df_all.Label == -1,
                       'original_target'] = df_all_sub['original_target'].tolist()
            
    
    
    df_all['original_target'] = df_all['original_target'].str.replace(
        "\\[|\\]|\\.|\\d+", "", regex=True)
    df_all = df_all.drop_duplicates(subset='original_target')
        

    df_all.drop(['original_target', 'sorted_Peptide', 'id', 'rank'], axis=1, inplace=True, errors='ignore')
    return(df_all)
#########################################################################################################


def custom_accuracy(y, y_pred, alpha=0.01, p=0.5, mult=1, BC1=1):
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
    qvals = uf.TDC_flex_c(decoy_wins, target_wins, BC1=BC1,
                          c=(mult*(1 - p) + 1)/(mult + 1), lam=(mult*(1 - p) + 1)/(mult + 1))  # doing TDC
    disc = sum((qvals[i] <= alpha) and (y[i] == 1)
               for i in range(len(y)))  # getting discoveries
    return(disc)
#########################################################################################################


def train_cv(labels, df, folds=3, Cs=[0.1, 1, 10], kernel='linear', degree=2, alpha=0.01, p=0.5, mult=1):
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

    if kernel not in ['linear', 'poly', 'rbf']:
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

    my_scorer = make_scorer(custom_accuracy, alpha=alpha, p=p, mult=mult,
                            greater_is_better=True, needs_threshold=True)

    if kernel == 'linear':
        grid = GridSearchCV(svm.LinearSVC(max_iter=int(1e7),dual=True), param_grid=param_grid,
                            cv=folds, scoring=my_scorer)
    else:
        grid = GridSearchCV(svm.SVC(), param_grid=param_grid,
                            cv=folds, scoring=my_scorer)
    grid.fit(df, labels)
    if max(grid.cv_results_['mean_test_score']) == 0:
        my_scorer = make_scorer(custom_accuracy, alpha=alpha, p=p, mult=mult, BC1=0,
                                greater_is_better=True, needs_threshold=True)

        if kernel == 'linear':
            grid = GridSearchCV(svm.LinearSVC(max_iter=int(1e7),dual=True), param_grid=param_grid,
                                cv=folds, scoring=my_scorer)
        else:
            grid = RandomizedSearchCV(svm.SVC(kernel=kernel), param_distributions=param_grid, n_iter=10,
                                      cv=folds, scoring=my_scorer)
        grid.fit(df, labels)
    
    while max(grid.cv_results_['mean_test_score']) == 0:
        alpha = alpha + 0.005
        my_scorer = make_scorer(custom_accuracy, alpha=alpha, BC1=0,
                                greater_is_better=True, needs_threshold=True)

        if kernel == 'linear':
            grid = GridSearchCV(svm.LinearSVC(max_iter=int(1e7),dual=True), param_grid=param_grid,
                                cv=folds, scoring=my_scorer)
        else:
            grid = RandomizedSearchCV(svm.SVC(kernel=kernel), param_distributions=param_grid, n_iter=10,
                                      cv=folds, scoring=my_scorer)

        grid.fit(df, labels)

    return(grid)
#########################################################################################################


def do_scale(df_all, df_extra=None):
    '''
    

    Parameters
    ----------
    df_all : Pandas dataframe
        Search file.
    df_extra : Pandas dataframe, optional
        Search file. The default is None.

    Returns
    -------
    Scales the features in df_all for SVM training. Applies the same transformation to df_extra.

    '''
    #scale non-binary features
    scale = StandardScaler()
    df_all.loc[:, ~(df_all.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins', 'trained']))] = scale.fit_transform(
        df_all.loc[:, ~(df_all.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins', 'trained']))])

    if df_extra is None:
        return(df_all, scale)
    else:
        df_extra.loc[:, ~(df_extra.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins', 'trained']))] = scale.transform(
            df_extra.loc[:, ~(df_extra.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins', 'trained']))])
        return(df_all, scale, df_extra)

################################################################################################


def do_svm(df_all, train_all, df_orig, folds=3, Cs=[0.1, 1, 10], total_iter=5, p=0.5, kernel='linear', alpha=0.01, train_alpha=0.01, degree=None, remove=None, top_positive=True, mult=1, initial_dir='XCorr'):
    '''
    

    Parameters
    ----------
    df_all : Pandas dataframe
        Search file with scaled features.
    train_all : Pandas dataframe
        Search file containing a subset of the Search file's decoy PSMs, scaled.
    df_orig : Pandas dataframe
        Search file with unscaled features.
    folds : Int, optional
        The number of folds for selection of class weights. The default is 3.
    Cs : List, optional
        Grid of class weights. The default is [0.01, 0.1, 1, 10].
    total_iter : Int, optional
        Number of SVM training rounds. The default is 10.
    p : Float, optional
        The fraction of decoys randomly sampled from df_all that are in train_all. If more than one decoy
        index used, this fraction is only taken with respect to the first decoy index. The default is 0.5.
    kernel : String, optional
        Type of kernel used for SVM. The default is 'linear'.
    alpha : Float, optional
        FDR threshold. The default is 0.01.
    train_alpha : Float, optional
        Training FDR threshold for selection of positive training set. The default is 0.01.
    degree : Int, optional
        Degree if polynomial kernel is used. The default is None.
    remove : List, optional
        A list of features to remove from training. The default is None.
    top_positive : Bool, optional
        Whether the top 1 PSM should be used for the positive set. The default is True.
    mult : Int, optional
         The number of decoy indices used. The default is 1.

    Returns
    -------
    Returns df_all, with the associated SVM scores, train_all with the associated SVM scores, and the learnt SVM model.

    '''
    train_alpha_init = train_alpha

    # ensure the indices of train_all belong to df_all
    train_targets = df_all[~(df_all.index.isin(train_all.index))].copy()
    train_targets.loc[:, 'Label'] = 1

    train_df = pd.concat([train_all, train_targets]).reset_index(drop=True)
    train_df = train_df.sample(frac=1).reset_index(drop=True)
    
    get_columns = train_df.columns[train_df.columns.str.contains(initial_dir, case=False)]
    
    if len(get_columns) > 0:
        initial_dir_case_insensitive = train_df.columns[train_df.columns.str.contains(initial_dir, case=False)][0]
    else:
        sys.exit("--initial_dir %s not detected. \n" %(initial_dir))
    
    train_df = train_df.sort_values(
        by=initial_dir_case_insensitive, ascending=False).reset_index(drop=True)
        
    #real df
    real_df = df_all[~(df_all.index.isin(train_all.index))
                     ].copy().reset_index(drop=True)
    df_orig = df_orig[~(df_orig.index.isin(train_all.index))
                      ].copy().reset_index(drop=True)
    
    #Preprocess dataframe
    SVM_train_labels = train_df['Label'].copy()
    
    SVM_train_features = train_df.drop(
        ['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins'], axis=1, errors='ignore').copy()
    
    if isinstance(remove, list):
        remove_list = [r for r in remove if r in SVM_train_features.columns]
        sys.stderr.write("Dropping the features from training: %s. \n" %(', '.join(remove_list)))
        logging.info("Dropping the features from training: %s." %(', '.join(remove_list)))        
        SVM_train_features.drop(remove_list, axis=1, inplace=True, errors='ignore')
    
    #Get rid of redundant features
    sds = SVM_train_features.apply(np.std, axis=0)
    SVM_train_features = SVM_train_features[SVM_train_features.columns[sds != 0]]

    columns_trained = SVM_train_features.columns

    #getting initial positive and negative set
    q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                           SVM_train_labels == 1, BC1=1, c=(mult*(1 - p) + 1)/(mult + 1), lam=(mult*(1 - p) + 1)/(mult + 1))
    
    if top_positive and ('rank' in SVM_train_features.columns):
        positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (
            SVM_train_features['rank'] == min(SVM_train_features['rank']))
    else:
        positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha)
    
    if sum(positive_set_indxs) == 0:
        logging.info("No initial positive set: removing the +1 penalty.")
        sys.stderr.write("No initial positive set: removing the +1 penalty.\n")
        q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                               SVM_train_labels == 1, BC1=0, c=(mult*(1 - p) + 1)/(mult + 1), lam=(mult*(1 - p) + 1)/(mult + 1))
        if top_positive and ('rank' in SVM_train_features.columns):
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (
                SVM_train_features['rank'] == min(SVM_train_features['rank']))
        else:
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha)
    
    while sum(positive_set_indxs) == 0:
        train_alpha = train_alpha + 0.005

        if top_positive and ('rank' in SVM_train_features.columns):
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (
                SVM_train_features['rank'] == min(SVM_train_features['rank']))
        else:
            positive_set_indxs = (SVM_train_labels == 1) & (
                q_vals <= train_alpha)
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
                        folds=folds, Cs=Cs, kernel=kernel, degree=degree, alpha=train_alpha, p=p, mult=mult)
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
                               SVM_train_labels == 1, BC1=1, c=(mult*(1 - p) + 1)/(mult + 1), lam=(mult*(1 - p) + 1)/(mult + 1))

        if top_positive and ('rank' in SVM_train_features.columns):
            positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(
                real_df.SpecId)) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
        else:
            positive_set_indxs = (SVM_train_labels == 1) & (
                q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId))
        negative_set_indxs = (SVM_train_labels == -1)
        
        if sum(positive_set_indxs) == 0:
            logging.info("No positive set: removing the +1 penalty.")
            sys.stderr.write("No positive set: removing the +1 penalty.\n")
            q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                                   SVM_train_labels == 1, BC1=0, c=(mult*(1 - p) + 1)/(mult + 1), lam=(mult*(1 - p) + 1)/(mult + 1))
            if top_positive and ('rank' in SVM_train_features.columns):
                positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(
                    real_df.SpecId)) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
            else:
                positive_set_indxs = (SVM_train_labels == 1) & (
                    q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId))

        while sum(positive_set_indxs) == 0:
            train_alpha = train_alpha + 0.005
            #getting initial positive and negative set
            q_vals = uf.TDC_flex_c(SVM_train_labels == -1,
                                   SVM_train_labels == 1, BC1=0, c=(mult*(1 - p) + 1)/(mult + 1), lam=(mult*(1 - p) + 1)/(mult + 1))
            if top_positive and ('rank' in SVM_train_features.columns):
                positive_set_indxs = (SVM_train_labels == 1) & (q_vals <= train_alpha) & (train_df.SpecId.isin(
                    real_df.SpecId)) & (SVM_train_features['rank'] == min(SVM_train_features['rank']))
            else:
                positive_set_indxs = (SVM_train_labels == 1) & (
                    q_vals <= train_alpha) & (train_df.SpecId.isin(real_df.SpecId))

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
        real_df_test = real_df.drop(
            ['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins'], axis=1, errors='ignore').copy()
        
        if isinstance(remove, list):
            remove_list = [r for r in remove if r in SVM_train_features.columns]     
            real_df_test.drop(remove_list, axis=1, inplace=True, errors='ignore')

        #Get rid of redundant features
        real_df_test = real_df_test[real_df_test.columns[real_df_test.columns.isin(
            SVM_train_features.columns)]]

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

        logging.info("Mean cross-validation power: %s." % (best_train_power))
        sys.stderr.write("Mean cross-validation power: %s.\n" % (best_train_power))

        logging.info("Std cross-validation power: %s." % (best_train_std))
        sys.stderr.write("Std cross-validation power: %s.\n" % (best_train_std))

    #using the last new_idx to report the discoveries
    df_orig['SVM_score'] = new_scores
    df_orig = df_orig.sort_values(by='SVM_score', ascending=False).reset_index(
        drop=True)  # sort by score

    q_val = uf.TDC_flex_c(
        df_orig.Label == -1, df_orig.Label == 1, c=1/(mult*(1 - p) + 1), lam=1/(mult*(1 - p) + 1))

    df_orig['q_val'] = q_val

    if isinstance(remove, list):
        train_all_test = train_all.drop(remove, axis=1, errors='ignore')
    else:
        train_all_test = train_all.copy()

    train_all_test = train_all_test.drop(
        ['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins'], axis=1, errors='ignore').copy()

    #Get rid of redundant features
    train_all_test = train_all_test[train_all_test.columns[train_all_test.columns.isin(
        SVM_train_features.columns)]]

    new_scores = grid.decision_function(train_all_test)

    train_all['SVM_score'] = new_scores

    return(train_power, train_std, true_power, df_orig, train_all, grid, columns_trained)


################################################################################################

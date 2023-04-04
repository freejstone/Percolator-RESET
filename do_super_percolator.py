#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 17:15:01 2023

@author: jackfreestone

This module performs the percolator algorithm with strict FDR control 
"""

import numpy as np
import pandas as pd


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


def peptide_level(narrow_file, open_file, peptide_list, score = 'TailorScore'):
    '''
    Parameters
    ----------
    narrow_file : string
        path to filtered narrow pin file.
    open_file : string
        path to filtered open pin file.
    peptide_list : string
        DESCRIPTION.
    score : TYPE, optional
        DESCRIPTION. The default is 'TailorScore'.

    Returns
    -------
    None.

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
    
    return(df_all)


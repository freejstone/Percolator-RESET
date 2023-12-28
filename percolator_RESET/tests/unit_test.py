#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 22:38:47 2023

@author: jackfreestone

Testing important utility/helper functions used in percolator RESET. Overall testing is accomplished by ensuring
that the percolator RESET completes successfully using a few example data sets via Github actions.

"""

import pandas as pd
from percolator_RESET import utility_functions as uf
from percolator_RESET import percolator_functions as pf

def test_TDC_flex_c():
    decoy_wins = pd.Series([False, False, True, True])
    decoy_wins1 = pd.Series([False, False, False, True])
    target_wins = pd.Series([True, True, False, False])
    assert uf.TDC_flex_c(decoy_wins, target_wins, BC1=1, c=1/2, lam=1/2).to_list() == [0.5, 0.5, 1.0, 1.0]
    assert uf.TDC_flex_c(decoy_wins1, target_wins, BC1=1, c=1/2, lam=1/2).to_list() == [0.5, 0.5, 0.5, 1.0]


def test_reverse_sequence():
    assert uf.reverse_sequence('ABCD') == 'CBAD'
    assert uf.reverse_sequence('ABCD[10.1]') == 'CBAD[10.1]'
    assert uf.reverse_sequence('A[9.2]BCD[10.1]') == 'CBA[9.2]D[10.1]'
    assert uf.reverse_sequence('n[8.3]A[9.2]BCD[10.1]') == 'C[8.3]BA[9.2]D[10.1]'


def test_check_n_term():
    sequences = pd.Series(['[1.2345]ABCDEFG', '[1.2345]BCDEFG'])
    assert uf.check_n_term(sequences).to_list() == ['A[1.2345]BCDEFG', 'B[1.2345]CDEFG']


def test_comet_pairing():
    search_file = './Percolator_RESET/tests/comet_open_top50.comet.pin'
    data_df = uf.read_pin(search_file)
    #removing flanking aa
    data_df['Peptide'] = data_df['Peptide'].str.extract(r'^[^.]*\.(.*?)\.[^.]*$', expand=False).fillna(data_df['Peptide'])
    data_df['Proteins'] = data_df['Proteins'].str.replace('\t', ',')
    data_dfs = [data_df]
    data_dfs, peptide_list_dfs = uf.comet_pairing(data_dfs, 'lnexpect')
    assert all(data_dfs[0].loc[data_dfs[0].Label == 1, 'Peptide'].isin(peptide_list_dfs[0].target))
    assert all(data_dfs[0].loc[data_dfs[0].Label == -1, 'Peptide'].isin(peptide_list_dfs[0].decoy))

def test_average_score():
    search_files = ['./Percolator_RESET/tests/narrow_50_sep.tide-search.target.pin', './Percolator_RESET/tests/narrow_51_sep.tide-search.target.pin']
    data_dfs = []
    for search_file in search_files:
        data_df = uf.read_pin(search_file)
        #removing flanking aa
        data_df['Peptide'] = data_df['Peptide'].str.extract(r'^[^.]*\.(.*?)\.[^.]*$', expand=False).fillna(data_df['Peptide'])
        data_df['Proteins'] = data_df['Proteins'].str.replace('\t', ',')
        data_dfs.append(data_df)
    average_scores = (data_dfs[0].TailorScore + data_dfs[1].TailorScore)/2
    data_dfs = uf.average_scores(data_dfs, ['TailorScore'])
    assert data_dfs[0].TailorScore.to_list() == average_scores.to_list()

def test_PSM_level():
    search_files = ['./Percolator_RESET/tests/narrow_50_sep.tide-search.target.pin', './Percolator_RESET/tests/narrow_50_sep.tide-search.decoy.pin']
    data_dfs = []
    for search_file in search_files:
        data_df = uf.read_pin(search_file)
        #removing flanking aa
        data_df['Peptide'] = data_df['Peptide'].str.extract(r'^[^.]*\.(.*?)\.[^.]*$', expand=False).fillna(data_df['Peptide'])
        data_df['Proteins'] = data_df['Proteins'].str.replace('\t', ',')
        data_dfs.append(data_df)
    results = pf.PSM_level(data_dfs[0], data_dfs[1], top=1)
    assert results.XCorr[0] == 0.0
    assert results.shape[0] == 1
    
def test_peptide_level():
    search_file = './Percolator_RESET/tests/comet_open_top50.comet2.pin'
    data_df = uf.read_pin(search_file)
    #removing flanking aa
    data_df['Peptide'] = data_df['Peptide'].str.extract(r'^[^.]*\.(.*?)\.[^.]*$', expand=False).fillna(data_df['Peptide'])
    data_df['Proteins'] = data_df['Proteins'].str.replace('\t', ',')
    data_dfs = [data_df]
    data_dfs, peptide_list_dfs = uf.comet_pairing(data_dfs, 'lnexpect')
    results = pf.peptide_level(data_dfs[0].copy(), peptide_list_dfs.copy(), True, 'lnExpect')
    assert results.shape[0] == 2
    
    search_file = './Percolator_RESET/tests/comet_open_top50.comet2.pin'
    data_df = uf.read_pin(search_file)
    #removing flanking aa
    data_df['Peptide'] = data_df['Peptide'].str.extract(r'^[^.]*\.(.*?)\.[^.]*$', expand=False).fillna(data_df['Peptide'])
    data_df['Proteins'] = data_df['Proteins'].str.replace('\t', ',')
    data_dfs = [data_df]
    data_dfs, peptide_list_dfs = uf.comet_pairing(data_dfs, 'lnexpect')
    data_dfs[0].loc[5, 'Peptide'] = peptide_list_dfs[0].loc[peptide_list_dfs[0].decoy == data_dfs[0].loc[0, 'Peptide'], 'target'].values[0]
    results = pf.peptide_level(data_dfs[0].copy(), peptide_list_dfs.copy(), True, 'lnExpect')
    assert results.shape[0] == 1


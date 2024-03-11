#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 15:31:40 2023

@author: jfre0619

Utility functions

"""
import os
import datetime
import platform
import numpy as np
import pandas as pd
import logging
import sys
import gzip
import re

#########################################################################################################


def print_info(command_line, output_dir, file_root, overwrite, search_files):
    #check if output directory exists, if not create and store log file there.
    ext = ".log.txt"
    level_message = '%(levelname)s: %(message)s'
    if os.path.isdir(output_dir):
        if os.path.exists(path=output_dir + "/" + file_root + ext) and overwrite:
            os.remove(output_dir + "/" + file_root + ext)
            logging.basicConfig(filename=output_dir + "/" + file_root + ext,
                                level=logging.DEBUG, format=level_message, force=True)
        elif os.path.exists(path=output_dir + "/" + file_root + ext) and not overwrite:
            log_file = output_dir + "/" + file_root + ext
            sys.exit("The file %s already exists and cannot be overwritten. Use --overwrite T to replace or choose a different output file name. \n" % (log_file))
        else:
            sys.stderr.write("check")
            logging.basicConfig(filename=output_dir + "/" + file_root + ext,
                                level=logging.DEBUG, format=level_message, force=True)
    else:
        sys.stderr.write("check1")
        os.mkdir(output_dir)
        logging.basicConfig(filename=output_dir + "/" + file_root + ext,
                            level=logging.DEBUG, format=level_message, force=True)

    #print CPU info
    logging.info('CPU: ' + str(platform.platform()))
    sys.stderr.write('CPU: ' + str(platform.platform()) + " \n")

    #print date time info
    logging.info(str(datetime.datetime.now()))
    sys.stderr.write(str(datetime.datetime.now()) + " \n")

    #print command used
    logging.info('Command used: ' + command_line)
    sys.stderr.write('Command used: ' + command_line + "\n")

    sys.stderr.write("Successfully read in arguments. \n")
    logging.info("Successfully read in arguments")

    if isinstance(search_files, list):
        for i in search_files:
            if not os.path.isfile(i):
                logging.info("%s does not exist." %(i))
                sys.exit("%s does not exist. \n" %(i)) 
    
#########################################################################################################


def TDC_flex_c(decoy_wins, target_wins, BC1=1, c=1/2, lam=1/2):
    '''
    Parameters
    ----------
    decoy_wins : boolean series
        indicates the position of the decoy wins.
    target_wins : boolean series
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
    fdps = np.minimum(1, ((BC1 + nDD) / np.maximum(nTD, 1)) * (c / (1-lam)))
    qvals = np.minimum.accumulate(fdps[::-1])[::-1]
    return(qvals)
#########################################################################################################


def read_pin(perc_file):
    """
    Read a Percolator tab-delimited file.

    Percolator input format (PIN) files and the Percolator result files
    are tab-delimited, but also have a tab-delimited protein list as the
    final column. This function parses the file and returns a DataFrame.

    Parameters
    ----------
    perc_file : str
        The file to parse.

    Returns
    -------
    pandas.DataFrame
        A DataFrame of the parsed data.
    """
    if str(perc_file).endswith(".gz"):
        fopen = gzip.open
    else:
        fopen = open

    with fopen(perc_file) as perc:
        cols = perc.readline().rstrip().split("\t")
        dir_line = perc.readline().rstrip().split("\t")[0]
        if dir_line.lower() != "defaultdirection":
            perc.seek(0)
            _ = perc.readline()

        psms = pd.concat((c for c in _parse_in_chunks(perc, cols)), copy=False)

    return psms
#########################################################################################################


def _parse_in_chunks(file_obj, columns, chunk_size=int(1e8)):
    """
    Parse a file in chunks

    Parameters
    ----------
    file_obj : file object
        The file to read lines from.
    columns : list of str
        The columns for each DataFrame.
    chunk_size : int
        The chunk size in bytes.

    Returns
    -------
    pandas.DataFrame
        The chunk of PSMs
    """
    constant_cols = ['SpecId', 'filename', 'Peptide', 'Proteins']
    
    while True:
        psms = file_obj.readlines(chunk_size)
        if not psms:
            break

        psms = [l.rstrip().split("\t", len(columns) - 1) for l in psms]
        psms = pd.DataFrame.from_records(psms, columns=columns)
        
        for column in columns:
            if column not in constant_cols:
                psms[column] = pd.to_numeric(psms[column]) #errors = 'ignore' is not allowed to be captured anymore!
        #psms = psms.apply(pd.to_numeric, errors="ignore")
        
        sys.stderr.write(psms.dtypes.to_string())
        
        yield psms
        
#########################################################################################################


def create_cluster(target_decoys, scale, original_discoveries, model, isolation_window, columns_trained, initial_dir):
    '''
    

    Parameters
    ----------
    target_decoys : Pandas dataframe
        top 1 PSMs, scaled.
    scale : StandardScaler
        For scaling and unscaling target_decoys
    original_discoveries : Pandas series
        Discovered peptides.
    model : SVM.SVC()
        SVM model.
    isolation_window : List
        Left/Right Isolation Window.
    columns_trained: List
        Names of the columns during SVM training
    Returns
    -------
    target PSMs associated to the discovered peptides that are the highest scoring in their associated mass bin.

    '''
    #get targets
    targets = target_decoys[target_decoys['Label'] == 1].copy()
    
    targets['original_target_sequence'] = targets['Peptide'].str.replace(
        "\\[|\\]|\\.|\\d+", "", regex=True)
        
    targets = targets[targets['original_target_sequence'].isin(
        original_discoveries)].copy()

    targets = targets.sample(frac=1)
    
    #getting charge
    targets['charge_temp'] = targets['SpecId'].apply(
        lambda x: int(x[-3]))
    
    #clustering the PSMs matched to the same peptide according to mass
    targets = targets.sort_values(
        by='ExpMass', ascending=True).reset_index(drop=True)
    targets["mass_plus"] = targets.groupby('original_target_sequence', group_keys=False).apply(
        lambda x: x.ExpMass.shift(1) + np.maximum(isolation_window[0]*x.charge_temp, isolation_window[1]*x.charge_temp.shift(1)))
    targets.loc[targets["mass_plus"].isna(), "mass_plus"] = -np.Inf
    targets["condition"] = targets["ExpMass"] > targets["mass_plus"]
    targets["cluster"] = targets.groupby(
        'original_target_sequence', group_keys=False).condition.cumsum()
    
    scaled_target = targets[targets.columns[~(targets.columns.isin(['SpecId', 'Label', 'filename', 'fileindx', 'ScanNr', 'Peptide', 'Proteins', 'trained', 'charge_temp', 'condition', 'mass_plus', 'cluster', 'original_target_sequence']))]].copy()
    
    scaled_target.loc[:,:] = scale.transform(scaled_target)
    
    new_scores = model.decision_function(scaled_target[scaled_target.columns[scaled_target.columns.isin(columns_trained)]])
    
    targets['SVM_score'] = new_scores
    
    targets = targets.sort_values(by=[initial_dir], ascending=False)
        
    #take best PSM according to cluster and sequence with modification
    targets = targets.drop_duplicates(subset=['Peptide', 'cluster'])
    
    targets = targets.drop(['charge_temp', 'condition', 'mass_plus', 'cluster', 'original_target_sequence'], axis = 1)

    return(targets)

#########################################################################################################


def score_PSMs(target_decoys, scale, model, columns_trained):
    '''

    Parameters
    ----------
    target_decoys : Pandas dataframe
        top 1 PSMs, scaled.
    scale : StandardScaler
        For scaling and unscaling target_decoys
    model : SVM.SVC()
        SVM model.
    columns_trained: List
        Names of the columns during SVM training
    Returns
    -------
    target PSMs are scored according to model.

    '''

    scaled_target_decoys = target_decoys[target_decoys.columns[~(target_decoys.columns.isin(['SpecId', 'Label', 'filename', 'fileindx', 'ScanNr', 'Peptide', 'Proteins']))]].copy()
    
    scaled_target_decoys.loc[:,:] = scale.transform(scaled_target_decoys)
    
    new_scores = model.decision_function(scaled_target_decoys[scaled_target_decoys.columns[scaled_target_decoys.columns.isin(columns_trained)]])
    
    target_decoys['SVM_score'] = new_scores

    return(target_decoys)
#########################################################################################################


def reverse_sequence(sequence):
    '''
    Parameters
    ----------
    sequence : Pandas Series
        List of sequence with variable modifications in tide-search format.

    Returns
    -------
    Recovers target sequences.
    '''
    results = re.findall(
        '[a-zA-Z]\\[\\d+.\\d+\\]\\[\\d+.\\d+\\]|[a-zA-Z]\\[\\d+.\\d+\\]|[a-zA-Z]\\[\\d+\\]|[a-zA-Z]', sequence)
    if 'n' in sequence:
        n_term_mass_mod_list = re.findall(r"(n\[\d+.\d+\])", sequence)
        if len(n_term_mass_mod_list) > 0:
            n_term_mass_mod = re.findall(r"(\d+.\d+)", n_term_mass_mod_list[0])
            results = results[1:]
            
    if 'c' in sequence:
        c_term_mass_mod_list = re.findall(r"(c\[\d+.\d+\])", sequence)
        if len(c_term_mass_mod_list) > 0:
            c_term_mass_mod = re.findall(r"(\d+.\d+)", c_term_mass_mod_list[0])
            results = results[:-1]


    results = results[-2::-1] + [results[-1]]

    if 'n' in sequence:
        results[0] = results[0] + '[' + n_term_mass_mod[0] + ']'

    if 'c' in sequence:
        results[-1] = results[-1] + '[' + c_term_mass_mod[0] + ']'

    results = ''.join(results)
    return(results)
#########################################################################################################


def check_n_term(sequences):
    '''
    Parameters
    ----------
    sequences : Pandas series
        List of sequence with variable modifications.

    Returns
    -------
    List of sequence with n-terminal modification properly relocated.

    '''

    n_term_bool = sequences.str.startswith('[')
    results = sequences[n_term_bool].str.split(pat=r'(\])', n=1, expand=True)

    results = results.apply(
        lambda x: x[2][0] + x[0][::] + x[1][::] + x[2][1:], axis=1)
    return(results)

#########################################################################################################


def comet_pairing(data_dfs, initial_dir):
    peptide_list_dfs = []
    for i in range(max(int(len(data_dfs)/2), 1)):
        
        j = int(i + len(data_dfs)/2)
        
        target_columns = data_dfs[i].Peptide[data_dfs[i].Label == 1].to_list()
        decoy_columns = data_dfs[j].Peptide[data_dfs[j].Label == -1].to_list()
        
        target_peptide_list = pd.Series(target_columns, name = 'target')
        decoy_peptide_list = pd.Series(decoy_columns, name = 'decoy')
        
        decoy_pairs = target_peptide_list.apply(reverse_sequence)
        target_pairs = decoy_peptide_list.apply(reverse_sequence)
        
        target_peptide_list = pd.DataFrame(zip(target_peptide_list, decoy_pairs), columns = ['target', 'decoy'])
        decoy_peptide_list = pd.DataFrame(zip(target_pairs, decoy_peptide_list), columns = ['target', 'decoy'])
        
        peptide_list_df = pd.concat([target_peptide_list, decoy_peptide_list])
        
        peptide_list_df = peptide_list_df.reset_index(drop = True)
        peptide_list_df['target'] = peptide_list_df['target'].str.replace("n", "").str.replace("c", "")
        peptide_list_df['decoy'] = peptide_list_df['decoy'].str.replace("n", "").str.replace("c", "")
        
        peptide_list_df.loc[peptide_list_df['target'].str.startswith('['), 'target'] = check_n_term(peptide_list_df['target'][peptide_list_df['target'].str.startswith('[')])
        peptide_list_df.loc[peptide_list_df['decoy'].str.startswith('['), 'decoy'] = check_n_term(peptide_list_df['decoy'][peptide_list_df['decoy'].str.startswith('[')])
        
        peptide_list_df.drop_duplicates(['target', 'decoy'], inplace = True)
        
        peptide_list_dfs.append(peptide_list_df)
        
        #to clean the peptide format so it is like tide/msfragger....where variable mods at terminals are not given
        
        data_dfs[i]['Peptide'] = data_dfs[i]['Peptide'].str.replace("n", "").str.replace("c", "")
        data_dfs[i].loc[data_dfs[i]['Peptide'].str.startswith('['), 'Peptide'] = check_n_term(data_dfs[i]['Peptide'][data_dfs[i]['Peptide'].str.startswith('[')])
        if (i != j):
            data_dfs[j]['Peptide'] = data_dfs[j]['Peptide'].str.replace("n", "").str.replace("c", "")
            data_dfs[j].loc[data_dfs[j]['Peptide'].str.startswith('['), 'Peptide'] = check_n_term(data_dfs[j]['Peptide'][data_dfs[j]['Peptide'].str.startswith('[')])
            
        if initial_dir.lower() == 'lnExpect'.lower():
            get_columns = data_dfs[i].columns[data_dfs[i].columns.str.contains(initial_dir, case=False)]
            
            if len(get_columns) > 0:
                initial_dir_case_insensitive = data_dfs[i].columns[data_dfs[i].columns.str.contains(initial_dir, case=False)][0]
            else:
                sys.exit("--initial_dir %s not detected. \n" %(initial_dir))
            
            
            sys.stderr.write("Taking the negative of lnExpect. \n")
            logging.info("Taking the negative of lnExpect.")
            data_dfs[i][initial_dir_case_insensitive] = -data_dfs[i][initial_dir_case_insensitive]
            if (i != j): 
                data_dfs[j][initial_dir_case_insensitive] = -data_dfs[j][initial_dir_case_insensitive]
    
    return(data_dfs, peptide_list_dfs)
#########################################################################################################


def average_scores(data_dfs, score):
    #averaging the two tailor PSMs.
    scores_list = [s for s in score if s in data_dfs[0].columns]
    
    if len(scores_list) > 0:
        data_dfs[0] = data_dfs[0].merge(
            data_dfs[1][['SpecId'] + scores_list], how='left', on='SpecId')
        
        sys.stderr.write("Averaging the scores: %s. \n" %(', '.join(scores_list)))
        logging.info("Averaging the scores: %s." %(', '.join(scores_list)))
        
        for s in scores_list:
            if s + '_x' in data_dfs[0].columns:             
                data_dfs[0][s] = (
                    data_dfs[0][s + '_x'] + data_dfs[0][s + '_y'])/2
                data_dfs[0].drop([s + '_x', s + '_y'],
                                 inplace=True, axis=1)
    data_dfs.pop(1)
    return(data_dfs)
#########################################################################################################
def report_psms(df_new, PSMs, get_psms, dynamic_competition, scale, model, isolation_window, columns_trained, initial_dir, output_dir, file_root):
    if df_new.shape[0] > 0 and get_psms and (dynamic_competition):
        sys.stderr.write(
            "Reporting all PSMs within each mass-cluster associated to a discovered peptide. \n")
        logging.info(
            "Reporting all PSMs within each mass-cluster associated to a discovered peptide. ")
        originally_discovered = df_new.loc[df_new.Label == 1, 'Peptide'].copy()
        originally_discovered = originally_discovered.str.replace(
            "\\[|\\]|\\.|\\d+", "", regex=True)
        df_extra = create_cluster(PSMs.copy(
        ), scale, originally_discovered, model, isolation_window, columns_trained, initial_dir)
        df_extra['originally_discovered'] = False
        df_new['originally_discovered'] = True
        
        #indicate whether the auxliary PSMs reported are above the associated SVM score threshold
        min_score = min(df_new.loc[df_new.Label == 1, 'SVM_score'])
        df_extra['above_threshold'] = False
        df_extra.loc[df_extra.SVM_score >= min_score, 'above_threshold'] = True
        df_new['above_threshold'] = True
    
        df_final = pd.concat([df_new, df_extra])
    
        # remove duplicate values from df_extra
        if any(df_final.columns.str.contains('filename')):
            df_final = df_final.drop_duplicates(['SpecId', 'filename', 'Peptide'])
        else:
            df_final = df_final.drop_duplicates(['SpecId', 'Peptide'])
        
        # write output
        df_final[df_final.Label == 1].to_csv(output_dir + "/" + file_root +
                                             ".psms.txt", header=True, index=False, sep='\t')

#########################################################################################################
    
def report_all_psms(PSMs, get_all_psms, scale, model, columns_trained, output_dir, file_root):
    if get_all_psms:
        PSMs = score_PSMs(PSMs.copy(), scale, model, columns_trained)
        PSMs.to_csv(output_dir + "/" + file_root +
                                             ".all_psms.txt", header=True, index=False, sep='\t')

#########################################################################################################

def report_decoys(df_new, train_all_unscale, train_all_new, report_decoys, output_dir, file_root):
    if report_decoys:
        decoys_final = df_new[df_new.Label == -1].reset_index(drop=True).copy()
        decoys_final['estimating_decoy'] = True
        decoys_final['training_decoy'] = False
        train_all_unscale['estimating_decoy'] = False
        train_all_unscale['training_decoy'] = True
        train_all_unscale['SVM_score'] = train_all_new.SVM_score
        decoys_final = pd.concat([decoys_final, train_all_unscale]).reset_index(drop = True).copy()
        decoys_final = decoys_final.sort_values(
            by='SVM_score', ascending=False).reset_index(drop=True)
        decoys_final.drop('q_val', axis=1, inplace=True, errors='ignore')
        # write output
        decoys_final.to_csv(output_dir + "/" + file_root +
                                         ".decoy_peptides.txt", header=True, index=False, sep='\t')
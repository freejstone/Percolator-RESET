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
#########################################################################################################


def print_info(command_line, output_dir, file_root, overwrite, search_files):
    #check if output directory exists, if not create and store log file there.
    if os.path.isdir(output_dir):
        if os.path.exists(path=output_dir + "/" + file_root + ".log.txt") and overwrite:
            os.remove(output_dir + "/" + file_root + ".log.txt")
            logging.basicConfig(filename=output_dir + "/" + file_root + ".log.txt",
                                level=logging.DEBUG, format='%(levelname)s: %(message)s', force=True)
        elif os.path.exists(path=output_dir + "/" + file_root + ".log.txt") and not overwrite:
            log_file = output_dir + "/" + file_root + ".log.txt"
            sys.exit("The file %s already exists and cannot be overwritten. Use --overwrite T to replace or choose a different output file name. \n" % (log_file))
        else:
            sys.stderr.write("check")
            logging.basicConfig(filename=output_dir + "/" + file_root + ".log.txt",
                                level=logging.DEBUG, format='%(levelname)s: %(message)s', force=True)
    else:
        sys.stderr.write("check1")
        os.mkdir(output_dir)
        logging.basicConfig(filename=output_dir + "/" + file_root + ".log.txt",
                            level=logging.DEBUG, format='%(levelname)s: %(message)s', force=True)

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

    if type(search_files) == list:
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
    while True:
        psms = file_obj.readlines(chunk_size)
        if not psms:
            break

        psms = [l.rstrip().split("\t", len(columns) - 1) for l in psms]
        psms = pd.DataFrame.from_records(psms, columns=columns)
        yield psms.apply(pd.to_numeric, errors="ignore")
#########################################################################################################


def create_cluster(target_decoys, scale, original_discoveries, model, isolation_window, columns_trained):
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
    
    targets['original_target_sequence'] = targets['Peptide']
    #targets['original_target_sequence'] = targets['Peptide'].str.replace(
    #    "\\[|\\]|\\.|\\d+", "", regex=True)
        
    targets = targets[targets['original_target_sequence'].isin(
        original_discoveries)].copy()

    targets = targets.sample(frac=1)
    
    #getting charge
    targets['charge'] = targets['SpecId'].apply(
        lambda x: int(x[-3]))
    
    #clustering the PSMs matched to the same peptide according to mass
    targets = targets.sort_values(
        by='ExpMass', ascending=True).reset_index(drop=True)
    targets["mass_plus"] = targets.groupby('original_target_sequence', group_keys=False).apply(
        lambda x: x.ExpMass.shift(1) + np.maximum(isolation_window[0]*x.charge, isolation_window[1]*x.charge.shift(1)))
    targets.loc[targets["mass_plus"].isna(), "mass_plus"] = -np.Inf
    targets["condition"] = targets["ExpMass"] > targets["mass_plus"]
    targets["cluster"] = targets.groupby(
        'original_target_sequence', group_keys=False).condition.cumsum()
    
    scaled_target = targets[targets.columns[~(targets.columns.isin(['SpecId', 'Label', 'filename', 'ScanNr', 'Peptide', 'Proteins', 'trained', 'charge', 'condition', 'mass_plus', 'cluster', 'original_target_sequence']))]].copy()
    
    scaled_target.loc[:,:] = scale.transform(scaled_target)
    
    new_scores = model.decision_function(scaled_target[scaled_target.columns[scaled_target.columns.isin(columns_trained)]])
    
    targets['SVM_score'] = new_scores
    
    if 'TailorScore' in targets.columns:
        targets = targets.sort_values(by=['TailorScore'], ascending=False)
    else:
        targets = targets.sort_values(by=['XCorr'], ascending=False)
        
    #take best PSM according to cluster and sequence with modification
    targets = targets.drop_duplicates(subset=['Peptide', 'cluster'])
    
    targets = targets.drop(['charge', 'condition', 'mass_plus', 'cluster', 'original_target_sequence'], axis = 1)

    return(targets)


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 13:40:07 2023

author: jackfreestone

This module performs the percolator RESET algorithm with strict FDR control using either a single or an extra decoy.
"""
import os
import time
import sys
import numpy as np
import pandas as pd
import random
import logging
from datetime import datetime
from percolator_RESET import utility_functions as uf
from percolator_RESET import percolator_functions as pf
from . import _version
__version__ = _version.get_versions()['version']

USAGE = """USAGE: python3 -m percolator_RESET [options] <search files> <target-decoy matchings>

  This script implements the percolator RESET algorithm with FDR control. The first input is a
  comma-separated list of the search file(s) in pin format to be used. The second input is a 
  comma-separated list of the associated target-decoy matches in the format of tide-index.
  
  For a single decoy database, provide <search files> with just the single concatenated search
  file (concat = T) and <target-decoy matchings> with the single target-decoy peptide list.
  
  In the case of using an extra decoy database, we need to provide 4 search files in <search files>:
  the two separate (concat = F) target search files and the two separate (concat = F) decoy
  search files. We also need to give the two target-decoy peptide lists in <target-decoy matchings>.
  Importantly, the two peptide lists need to be in order of the decoys provided. So if "_1" represents
  the first decoy database, and "_2" represents the second decoy database then:
      
  <search files> = 'target_1,target_2,decoy_1,decoy_2'
  <target-decoy matchings> = 'peptide_list_1,peptide_list_2'
  
  As an example:

      python3 -m percolator_RESET search_file1.txt peptide_list1.txt
      
      python3 -m percolator_RESET target_file1.txt,target_file2.txt,decoy_file1.txt,decoy_file2.txt peptide_list_1.txt,peptide_list_2.txt
  
  -----------------------------------------------------------------------------------------------------------------------------------------------    
  
  It is possible to run percolator RESET without explicitly pairing the targets and decoys together, thus alleviating the need for a peptide list.
  
  For the single decoy database version of RESET:
      
      python3 -m percolator_RESET --pair F search_file1.txt
      
  For the two decoy database version of RESET:
      
      python3 -m percolator_RESET --pair F --mult 2 search_file1.txt
      
  The --mult option here is important since the module will infer that you only used one decoy database.
  
  -----------------------------------------------------------------------------------------------------------------------------------------------
  
  Lastly, you may want to run percolator RESET having already completed the competition step externally. In this case for single decoy RESET:
      
      python3 -m percolator_RESET --dynamic_competition F search_file1.txt
      
  and for two decoy RESET:
      
      python3 -m percolator_RESET --dynamic_competition F --mult 2 search_file1.txt
            
  Options:
      
    --FDR_threshold <float> The FDR threshold. Default = 0.01.
    --folds <int> The number of folds for determining the class weights. Default = 3.
    --total_iter <int> The number of SVM training iterations. Default = 5.
    --train_FDR_threshold <float> The FDR threshold used to define positive training set. Default = 0.01.
    --output_dir <str> The path to the output directory. Default = '.'
    --report_decoys <T/F> A boolean indicating whether you would like the estimating decoys to be printed as well. Default = T.
    --file_root <str> The name of the file root when outputting results. Default = 'FDR_percolator'.
    --remove <str> A comma-separated list of features to remove prior to SVM trainings. Default = 'enzInt'.
    --overwrite <T/F> A boolean determining whether do_FDR_percolator.py should overwrite the present files in the directory. Default = F.
    --seed <int> Random seed. Default = int(datetime.now().timestamp()).
    --p_init <float> The proportion of decoys to use as the informative decoy set. Default = 0.5
    --get_psms <T/F> Prints out the relevant psms with distinct delta masses/variable mods associated with the discovered peptides. Default = F.
    --isolation_window <str> A comma-separated pair of numbers that describe the lower and upper isolation window. used for when get_psms = T. Default = 2,2
    --pair <T/F> A boolean determining whether explicit head-to-head competition using target-decoy pairing should be used. Default = T.
    --reverse <T/F> For Comet, a boolean indicating whether decoy peptides are obtained by reversing the target peptides (keeping C-teriminal fixed). If T, peptide pairing can be done without a user supplied peptide list. Default = F.
    --initial_dir <str> A string indicating which column will be the initial direction that defines the positive set for training and the peptide-level competition. Default = 'XCorr'.
    --score <str> A comma-separated list of calibrated scores that need averaging for the two-decoy database version of RESET (the same target PSM may have different scores depending on the decoy database used). Default = ''.
    --dynamic_competition <T/F> A boolean determining whether dynamic competition is first required. If not, it is assumed that the competition has been completed externally. If False, only a single search file in pin format is required. Default = T.
    --mult <int> The number of decoy databases used (this is usually inferred from the number of search files given but can be overidden if a single search file is provided but mult = 2 is required). Default = None.
"""


def main():
    global USAGE

    start_time = time.time()

    # Set default values for parameters.
    FDR_threshold = 0.01
    folds = 3
    total_iter = 5
    train_FDR_threshold = 0.01
    output_dir = '.'
    report_decoys = False
    file_root = 'FDR_percolator'
    remove = ['enzInt','ExpMass', 'CalcMass']
    overwrite = False
    # importantly the seed needs be random and not set at some value
    seed = int(datetime.now().timestamp())
    p_init_default = True
    get_psms = False
    get_all_psms = False
    isolation_window = [2, 2]
    pair = True
    initial_dir = 'XCorr'
    score = ['']
    dynamic_competition = True
    mult = None
    reverse= False
    command_line = ' '.join(sys.argv)

    # Parse the command line.
    sys.argv = sys.argv[1:]
    while (any('--' in string for string in sys.argv)):
        next_arg = sys.argv[0]
        sys.argv = sys.argv[1:]
        if (next_arg == "--FDR_threshold"):
            FDR_threshold = float(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--folds"):
            folds = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--total_iter"):
            total_iter = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--train_FDR_threshold"):
            train_FDR_threshold = float(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--output_dir"):
            output_dir = str(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--report_decoys"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                report_decoys = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                report_decoys = False
            else:
                sys.stderr.write("Invalid argument for --report_decoys")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--file_root"):
            file_root = str(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--remove"):
            remove = str(sys.argv[0]).split(',')
            sys.argv = sys.argv[1:]
        elif (next_arg == "--overwrite"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                overwrite = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                overwrite = False
            else:
                sys.stderr.write("Invalid argument for --overwrite")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == '--seed'):
            seed = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == '--p_init'):
            p_init_default = False
            p_init = float(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--get_psms"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                get_psms = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                get_psms = False
            else:
                sys.stderr.write("Invalid argument for --get_psms")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--get_all_psms"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                get_all_psms = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                get_all_psms = False
            else:
                sys.stderr.write("Invalid argument for --get_all_psms")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--isolation_window"):
            isolation_window = str(sys.argv[0]).split(',')
            isolation_window = [float(c) for c in isolation_window]
            sys.argv = sys.argv[1:]
        elif (next_arg == "--pair"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                pair = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                pair = False
            else:
                sys.stderr.write("Invalid argument for --pair")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--initial_dir"):
            initial_dir = str(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--score"):
            score = str(sys.argv[0]).split(',')
            sys.argv = sys.argv[1:]
        elif (next_arg == "--dynamic_competition"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                dynamic_competition = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                dynamic_competition = False
            else:
                sys.stderr.write("Invalid argument for --dynamic_competition")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--reverse"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                reverse = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                reverse = False
            else:
                sys.stderr.write("Invalid argument for --reverse")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == '--mult'):
            mult = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        else:
            sys.stderr.write("Invalid option (%s)" % next_arg)
            sys.exit(1)
    if (len(sys.argv) == 2):
        search_files = str(sys.argv[0]).split(',')
        td_lists = str(sys.argv[1]).split(',')
    elif ((len(sys.argv) == 1 and not dynamic_competition) or (len(sys.argv) == 1 and not pair) or (len(sys.argv) == 1 and reverse)):
        search_files = str(sys.argv[0]).split(',')
    else:
        sys.stderr.write('Version: ' + str(__version__) + " \n")
        sys.stderr.write(USAGE)
        sys.exit(1)
        
    ######### CODE STARTS HERE #########
    
    #setting seed for reproducibility
    if isinstance(seed, int):
        random.seed(seed)
        np.random.seed(seed)

    #print version
    logging.info('Version: ' + str(__version__))
    sys.stderr.write('Version: ' + str(__version__) + " \n")

    #print meta information, checking directory and printing warnings
    uf.print_info(command_line, output_dir, file_root, overwrite, search_files)

    sys.stderr.write("Using seed: %s. \n" % seed)
    logging.info("Using seed: %s." % seed)

    single_decoy = (len(search_files) == 1)
    
    if not dynamic_competition:
        pair = False
    
    sys.stderr.write("Reading in search file(s). \n")
    logging.info("Reading in search file(s).")
    
    #read in search files
    data_dfs = []
    for search_file in search_files:
        data_df = uf.read_pin(search_file)
        #removing flanking aa
        data_df['Peptide'] = data_df['Peptide'].str.extract(r'^[^.]*\.(.*?)\.[^.]*$', expand=False).fillna(data_df['Peptide'])
        data_df['Proteins'] = data_df['Proteins'].str.replace('\t', ',')
        data_dfs.append(data_df)
        
    if not dynamic_competition:
        sys.stderr.write("Skipping dynamic level competition (and assuming competition has already been performed). \n")
        logging.info("Skipping dynamic level competition (and assuming competition has already been performed).")

    if pair and not reverse:
        sys.stderr.write("Reading in peptide list(s). \n")
        logging.info("Reading in search file(s).")
        peptide_list_dfs = []
        for td_list in td_lists:
            peptide_list_df = pd.read_table(td_list)
            if 'decoy(s)' in peptide_list_df.columns:
                peptide_list_df.rename(columns={'decoy(s)': 'decoy'}, inplace=True)
            peptide_list_df.drop_duplicates(inplace=True)
            peptide_list_dfs.append(peptide_list_df)
    elif pair and reverse:
        sys.stderr.write("Pairing target and decoy peptides. \n")
        logging.info("Pairing target and decoy peptides.")
        #do pairing
        data_dfs, peptide_list_dfs = uf.comet_pairing(data_dfs, initial_dir)
    else:
        if dynamic_competition:
            sys.stderr.write("No peptide-list pairing provided. Not performing explicit head-to-head competition. \n")
            logging.info("No peptide-list pairing provided. Not performing explicit head-to-head competition.")
        else:
            df_all = data_dfs[0]
            
    if mult is None and single_decoy:
        mult = 1
    elif mult is None and not single_decoy:
        mult = 2
    
    if single_decoy:
        data_df = data_dfs[0]
    if not single_decoy:
        data_dfs = uf.average_scores(data_dfs, score)
        #do PSM level competition involving all PSMs
        data_df = pf.PSM_level(data_dfs[0].copy(), data_dfs[1:].copy(), top=1)
    
    if pair:
        #removing low complexity targets that do not produce a decoy
        low_complex_targets = peptide_list_dfs[0].target[peptide_list_dfs[0].decoy.isna(
        )]
        data_df = data_df[~data_df.Peptide.isin(
            low_complex_targets)].reset_index(drop=True)
        
    if dynamic_competition:
        #doing multi peptide level competition
        if pair:
            df_all = pf.peptide_level(
                data_df.copy(), peptide_list_dfs.copy(), pair, initial_dir)
        else:
            df_all = pf.peptide_level(
                data_df.copy(), None, pair, initial_dir)
    
    PSMs = data_df.copy()
    PSMs['rank'] = PSMs['SpecId'].apply(
        lambda x: int(x[-1]))
    PSMs = PSMs[PSMs['rank'] == 1].reset_index(drop=True)
    PSMs.drop('rank', inplace = True, axis = 1)
    
    #applying scaling
    df_all_scale, scale = pf.do_scale(df_all.copy())

    if p_init_default:
        p_init = 0.5

    #create target-decoys at pseudolevel
    rand_indxs = np.random.choice([True, False], replace = True, size = sum((df_all_scale['Label']
                                  == -1)), p = [p_init, 1 - p_init])
    
    train_all = df_all_scale.loc[(df_all_scale['Label']
                                  == -1)].copy()
    train_all = train_all.loc[rand_indxs].copy()
    
    train_all_unscale = df_all.loc[(df_all['Label']
                                  == -1)].copy()
    train_all_unscale = train_all_unscale.loc[rand_indxs].copy()

    #do SVM
    df_new, train_all_new, model, columns_trained = pf.do_svm(df_all_scale.copy(), train_all.copy(), df_all.copy(), folds=folds,
            p=p_init, total_iter=total_iter, alpha=FDR_threshold, train_alpha=train_FDR_threshold, remove = remove, mult=mult, initial_dir=initial_dir)

    df_new = df_new.loc[(df_new.q_val <= FDR_threshold) | (df_new.Label == -1)]
    
    coefficients = np.concatenate((model.best_estimator_.coef_[0], model.best_estimator_.intercept_))
    coefficients = pd.Series(coefficients)
    coefficients.index = columns_trained.append(pd.Index(['intercept']))
    
    logging.info('Final SVM feature coefficients:')
    sys.stderr.write('Final SVM feature coefficients: \n')
    logging.info(coefficients.to_string())
    sys.stderr.write(coefficients.to_string() + "\n")
    logging.info('Features that are constant are dropped.')
    sys.stderr.write('Features that are constant are dropped. \n')
    
    #write output folder
    if output_dir != '.':
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
    
    #report PSMs using clustering algorithm
    uf.report_psms(df_new, PSMs, get_psms, dynamic_competition, scale, model, isolation_window, columns_trained, initial_dir, output_dir, file_root)
    #report all PSMs
    uf.report_all_psms(PSMs, get_all_psms, scale, model, columns_trained, output_dir, file_root)
    #report discovered peptides
    df_new[df_new.Label == 1].to_csv(output_dir + "/" + file_root +
                                         ".peptides.txt", header=True, index=False, sep='\t')
    #report decoys
    uf.report_decoys(df_new, train_all_unscale, train_all_new, report_decoys, output_dir, file_root)
    
    end_time = time.time()

    logging.info("Elapsed time: " +
                 str(round(end_time - start_time, 2)) + " s")
    sys.stderr.write("Elapsed time: " +
                     str(round(end_time - start_time, 2)) + " s \n")


if __name__ == "__main__":
    main()

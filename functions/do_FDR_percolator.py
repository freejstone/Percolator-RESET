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
import utility_functions as uf
import percolator_functions as pf

USAGE = """USAGE: python3 do_FDR_percolator.py [options] <search files> <target-decoy matchings>

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

      python3 do_FDR_percolator.py search_file1.txt peptide_list1.txt
      
      python3 do_FDR_percolator.py target_file1.txt,target_file2.txt,decoy_file1.txt,decoy_file2.txt peptide_list_1.txt,peptide_list_2.txt
  
  -----------------------------------------------------------------------------------------------------------------------------------------------    
  
  It is possible to run percolator RESET without explicity pairing the targets and decoys together, thus allievating the need for a peptide list.
  
  For the single decoy database version of RESET:
      
      python3 do_FDR_percolator.py --pair F search_file1.txt
      
  For the two decoy database version of RESET:
      
      python3 do_FDR_percolator.py --pair F --mult 2 search_file1.txt
      
  The --mult option here is important since the module will infer that you only used one decoy database.
  
  -----------------------------------------------------------------------------------------------------------------------------------------------
  
  Lastly, you may want to run percolator RESET having already completed the competition step externally. In this case for single decoy RESET:
      
      python3 do_FDR_percolator.py --dynamic_competition F search_file1.txt
      
  and for two decoy RESET:
      
      python3 do_FDR_percolator.py --dynamic_competition F --mult 2 search_file1.txt
            
  Options:
      
    --FDR_threshold <float> The FDR threshold. Default = 0.01.
    --narrow <T/F> A boolean determining whether the search files used a narrow or open precursor tolerance. Default = T. 
    --folds <int> The number of folds for determining the class weights. Default = 3.
    --total_iter <int> The number of SVM training iterations. Default = 5.
    --train_FDR_threshold <float> The FDR threshold used to define positive training set. Default = 0.01.
    --output_dir <str> The path to the output directory. Default = '.'
    --report_decoys <T/F> A boolean indicating whether you would like the estimating decoys to be printed as well. Default = T.
    --file_root <str> The name of the file root when outputting results. Default = 'FDR_percoaltor'.
    --remove <str> A comma-separated list of features to remove prior to SVM trainings. Default = 'enzInt'.
    --overwrite <T/F> A boolean determining whether do_FDR_percolator.py should ovewrite the present files in the directory. Default = F.
    --seed <int> Random seed. Default = int(datetime.now().timestamp()).
    --p_init <float> The proportion of decoys to use as the informative decoy set. Default = 0.5
    --get_psms <T/F> Prints out the relevant psms with distinct delta masses/variable mods associated with the discovered peptides. Only relevant for --narrow F. Default = F.
    --isolation_window <str> A comma-separated pair of numbers that describe the lower and upper isolation window. used for when get_psms = T. Default = 2,2
    --pair <T/F> A boolean determining whether explicit head-to-head competition using target-decoy pairing should be used. Default = T.
    --initial_dir <str> A string indicating which column will be the initial direction that defines the positive set for training and the peptide-level competition. Default = 'XCorr'.
    --scores <str> A comma-separated list of calibrated scores that need averaging for the two-decoy database version of RESET (the same target PSM may have different scores depending on the decoy database used). Default = ''.
    --dynamic_competition <T/F> A boolean determining whether dynamic competition is first required. If not, it is assumed that the competition has been completed externally. If False, only a single search file in pin format is required. Default = T.
    --mult <int> The number of decoy databases used (this is usually inferred from the number of search files given but can be overidden if a single search file is provided but mult = 2 is required). Default = None.
"""


def main():
    global USAGE

    start_time = time.time()

    # Set default values for parameters.
    FDR_threshold = 0.01
    narrow = True
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
    scores = ['']
    dynamic_competition = True
    mult = None
    command_line = ' '.join(sys.argv)

    # Parse the command line.
    sys.argv = sys.argv[1:]
    while (any('--' in string for string in sys.argv)):
        next_arg = sys.argv[0]
        sys.argv = sys.argv[1:]
        if (next_arg == "--FDR_threshold"):
            FDR_threshold = float(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--narrow"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                narrow = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                narrow = False
            else:
                sys.stderr.write("Invalid argument for --narrow")
                sys.exit(1)
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
        elif (next_arg == "--scores"):
            scores = str(sys.argv[0]).spli(',')
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
        elif (next_arg == '--mult'):
            mult = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        else:
            sys.stderr.write("Invalid option (%s)" % next_arg)
            sys.exit(1)
    if (len(sys.argv) == 2):
        search_files = str(sys.argv[0]).split(',')
        td_lists = str(sys.argv[1]).split(',')
    elif ((len(sys.argv) == 1 and not dynamic_competition) or (len(sys.argv) == 1 and not pair)):
        search_files = str(sys.argv[0]).split(',')
    else:
        #sys.stderr.write('Version: ' + str(__version__) + " \n")
        sys.stderr.write(USAGE)
        sys.exit(1)

    #setting seed for reproducibility
    if type(seed) == int:
        random.seed(seed)
        np.random.seed(seed)

    #print meta information, checking directory and printing warnings
    uf.print_info(command_line, output_dir, file_root, overwrite, search_files)

    sys.stderr.write("Using seed: %s. \n" % seed)
    logging.info("Using seed: %s." % seed)

    single_decoy = (len(search_files) == 1)
    
    sys.stderr.write("Reading in search file(s). \n")
    logging.info("Reading in search file(s).")
    
    data_dfs = []
    for search_file in search_files:
        data_df = uf.read_pin(search_file)
        #removing flanking aa
        data_df['Peptide'] = data_df['Peptide'].str.extract(r'^[^.]*\.(.*?)\.[^.]*$', expand=False).fillna(data_df['Peptide'])
        data_dfs.append(data_df)
        
    if not dynamic_competition:
        sys.stderr.write("Skipping dynamic level competition (and assuming competition has already been performed). \n")
        logging.info("Skipping dynamic level competition (and assuming competition has already been performed).")

    if pair:
        sys.stderr.write("Reading in peptide list(s). \n")
        logging.info("Reading in search file(s).")
        peptide_list_dfs = []
        for td_list in td_lists:
            peptide_list_df = pd.read_table(td_list)
            if 'decoy(s)' in peptide_list_df.columns:
                peptide_list_df.rename(columns={'decoy(s)': 'decoy'}, inplace=True)
            peptide_list_df.drop_duplicates(inplace=True)
            peptide_list_dfs.append(peptide_list_df)
            
    else:
        if dynamic_competition:
            sys.stderr.write("No peptide-list pairing provided. Not performing explicit head-to-head competition. \n")
            logging.info("No peptide-list pairing provided. Not performing explicit head-to-head competition.")
        df_all = data_dfs[0]

    if single_decoy:
        
        if mult == None:
            mult = 1
        
        if pair:
            #removing low complexity targets that do not produce a decoy
            low_complex_targets = peptide_list_dfs[0].target[peptide_list_dfs[0].decoy.isna(
            )]
            data_dfs[0] = data_dfs[0][~data_dfs[0].Peptide.isin(
                low_complex_targets)].reset_index(drop=True)

        if dynamic_competition:
            #doing peptide level competition
            if pair:
                df_all = pf.peptide_level(
                    data_dfs[0].copy(), peptide_list_dfs[0].copy(), narrow, pair, initial_dir)
            else:
                df_all = pf.peptide_level(
                    data_dfs[0].copy(), None, narrow, pair, initial_dir)
        
        
        PSMs = data_dfs[0].copy()
        PSMs['rank'] = PSMs['SpecId'].apply(
            lambda x: int(x[-1]))
        PSMs = PSMs[PSMs['rank'] == 1].reset_index(drop=True)


    else:
        
        if mult == None:
            mult = 2
        
        #averaging the two tailor PSMs.
        scores_list = [s for s in scores if s in data_dfs[0].columns]
        
        data_dfs[0] = data_dfs[0].merge(
            data_dfs[1][['SpecId'] + scores_list], how='left', on='SpecId')
        
        sys.stderr.write("Averaging the scores: %s. \n" %(', '.join(scores_list)))
        logging.info("Averaging the scores." %(', '.join(scores_list)))
        
        for s in scores_list:
            if s + '_x' in data_dfs[0].columns:             
                data_dfs[0][s] = (
                    data_dfs[0][s + '_x'] + data_dfs[0][s + '_y'])/2
                data_dfs[0].drop([s + '_x', s + '_y'],
                                 inplace=True, axis=1)
        data_dfs.pop(1)
        
        #do PSM level competition involving all PSMs
        data_df = pf.PSM_level(data_dfs[0].copy(), data_dfs[1:].copy(), top=1)
        
        if pair:
            #removing low complexity targets that do not produce a decoy
            low_complex_targets = peptide_list_dfs[0].target[peptide_list_dfs[0].decoy.isna(
            )]
            data_df = data_df[~data_df.Peptide.isin(
                low_complex_targets)].reset_index(drop=True)

        PSMs = data_df.copy()
        PSMs['rank'] = PSMs['SpecId'].apply(
            lambda x: int(x[-1]))
        PSMs = PSMs[PSMs['rank'] == 1].reset_index(drop=True)
        
        if dynamic_competition:
            #doing multi peptide level competition
            if pair:
                df_all = pf.peptide_level(
                    data_df.copy(), peptide_list_dfs.copy(), narrow, pair, initial_dir)
            else:
                df_all = pf.peptide_level(
                    data_df.copy(), None, narrow, pair, initial_dir)

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
    train_power, std_power, true_power, df_new, train_all_new, model, columns_trained = pf.do_svm(df_all_scale.copy(), train_all.copy(), df_all.copy(), folds=folds, Cs=[
        0.1, 1, 10], p=p_init, total_iter=total_iter, alpha=FDR_threshold, train_alpha=train_FDR_threshold, remove = remove, mult=mult, initial_dir=initial_dir)

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
    
    if df_new.shape[0] > 0 and get_psms and (narrow == False) and (dynamic_competition):
        sys.stderr.write(
            "Reporting all PSMs within each mass-cluster associated to a discovered peptide. \n")
        logging.info(
            "Reporting all PSMs within each mass-cluster associated to a discovered peptide. ")
        originally_discovered = df_new.loc[df_new.Label == 1, 'Peptide'].copy()
        originally_discovered = originally_discovered.str.replace(
            "\\[|\\]|\\.|\\d+", "", regex=True)
        df_extra = uf.create_cluster(PSMs.copy(
        ), scale, originally_discovered, model, isolation_window, columns_trained, initial_dir)
        df_extra['originally_discovered'] = False
        df_new['originally_discovered'] = True
        
        #indicate whether the auxliary PSMs reported are above the associated SVM score threshold
        min_score = min(df_new.SVM_score)
        df_extra['above_threshold'] = False
        df_extra.loc[df_extra.SVM_score >= min_score, 'above_threshold'] = True
        df_new['above_threshold'] = True

        df_final = pd.concat([df_new, df_extra])

        # remove duplicate values from df_extra
        df_final = df_final.drop_duplicates(['SpecId', 'filename', 'Peptide'])
        
        # write output
        df_final[df_final.Label == 1].to_csv(output_dir + "/" + file_root +
                                             ".psms.txt", header=True, index=False, sep='\t')
        
    if get_all_psms:
        PSMs = uf.score_PSMs(PSMs.copy(), scale, model, columns_trained)
        PSMs.to_csv(output_dir + "/" + file_root +
                                             ".all_psms.txt", header=True, index=False, sep='\t')
    
    # write output
    df_new[df_new.Label == 1].to_csv(output_dir + "/" + file_root +
                                         ".peptides.txt", header=True, index=False, sep='\t')
        
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
        
    end_time = time.time()

    logging.info("Elapsed time: " +
                 str(round(end_time - start_time, 2)) + " s")
    sys.stderr.write("Elapsed time: " +
                     str(round(end_time - start_time, 2)) + " s \n")


if __name__ == "__main__":
    main()

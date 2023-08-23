#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 13:40:07 2023

author: jackfreestone

This module performs the percolator algorithm with strict FDR control using either a single or an extra decoy.
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

  This script implements the percolator algorithm with FDR control. The first input is a
  comma-separated list of the search file(s) in pin format to be used. The second input is a 
  comma-separated list of the associated target-decoy matches produced by tide-index.
  
  For a single decoy database, provide <search files> with just the single concatenated search
  file (concat = T) and <target-decoy matchings> with the single target-decoy peptide list.
  
  In the case of using an extra decoy database, we need to provide 4 search files in <search files>:
  the two separate (concat = F) target search files and the two separate (concat = F) decoy
  search files. We also need to give the two target-decoy peptide lists in <target-decoy matchings>.
  Importantly, the two peptide lists need to be in order of the decoys provided. So if "_1" represents
  the first decoy database, and "_2" represents the second decoy database then:
      
  <search files> = 'target_1,target_2,decoy_1,decoy_2'
  <target-decoy matchings> = 'peptide_list_1,peptide_list_2'
            
  Options:
      
    --FDR_threshold <float> The FDR threshold. Default = 0.01.
    --folds <int> The number of folds for determining the class weights. Default = 3.
    --total_iter <int> The number of SVM training iterations. Default = 5.
    --train_FDR_threshold <float> The FDR threshold used to define positive training set. Default = 0.01.
    --output_dir <str> The path to the output directory. Default = '.'
    --file_root <str> The name of the file root when outputting results. Default = 'FDR_percoaltor'.
    --remove <str> A comma-separated list of features to remove prior to SVM trainings. Default = 'enzInt'.
    --overwrite <T/F> A boolean determining whether do_FDR_percolator.py should ovewrite the present files in the directory. Default = F.
    --seed <int> Random seed. Default = int(datetime.now().timestamp()).
    --p_init <float> The proportion of decoys to use as the informative decoy set. Default = 0.5
    --get_psms <T/F> Prints out the relevant psms with distinct delta masses/variable mods associated with the discovered peptides. Default = F.
    --isolation_window <str> A comma-separated pair of numbers that describe the lower and upper isolation window. used for when get_psms = T. Default = 2,2
    
    
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
    file_root = 'FDR_percolator'
    remove = ['enzInt']
    peptide_level = False
    overwrite = False
    # importantly the seed needs be random and not set at some value
    seed = int(datetime.now().timestamp())
    p_init_default = True
    get_psms = False
    isolation_window = [2, 2]
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
        elif (next_arg == "--file_root"):
            file_root = str(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--remove"):
            remove = str(sys.argv[0]).split(',')
            sys.argv = sys.argv[1:]
        elif (next_arg == "--peptide_level"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                peptide_level = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                peptide_level = False
            else:
                sys.stderr.write("Invalid argument for --overwrite")
                sys.exit(1)
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
        elif (next_arg == "--isolation_window"):
            isolation_window = str(sys.argv[0]).split(',')
            isolation_window = [float(c) for c in isolation_window]
            sys.argv = sys.argv[1:]
        else:
            sys.stderr.write("Invalid option (%s)" % next_arg)
            sys.exit(1)
    if (len(sys.argv) == 2):
        search_files = str(sys.argv[0]).split(',')
        td_lists = str(sys.argv[1]).split(',')
    else:
        #sys.stderr.write('Version: ' + str(__version__) + " \n")
        sys.stderr.write(USAGE)
        sys.exit(1)

    #setting seed for reproducibility
    if type(seed) == int:
        random.seed(seed)
        np.random.seed(seed)

        sys.stderr.write("Using seed: %s. \n" % seed)
        logging.info("Using seed: %s." % seed)

    #print meta information, checking directory and printing warnings
    uf.print_info(command_line, output_dir, file_root, overwrite, search_files)

    

    single_decoy = (len(search_files) == 1)

    if single_decoy:
        
        sys.stderr.write("Reading in search file(s) and peptide list(s). \n")
        logging.info("Reading in search file(s) and peptide list(s).")
        
        data_dfs = []
        peptide_list_dfs = []
        for search_file in search_files:
            data_df = uf.read_pin(search_file)
            #removing flanking aa
            data_df['Peptide'] = data_df['Peptide'].apply(
                lambda x: x[2:(len(x) - 2)])
            data_dfs.append(data_df)
        for td_list in td_lists:
            peptide_list_df = pd.read_table(td_list)
            if 'decoy(s)' in peptide_list_df.columns:
                peptide_list_df.rename(columns={'decoy(s)': 'decoy'}, inplace=True)
            peptide_list_df.drop_duplicates(inplace=True)
            peptide_list_dfs.append(peptide_list_df)
        
        #removing low complexity targets that do not produce a decoy
        low_complex_targets = peptide_list_dfs[0].target[peptide_list_dfs[0].decoy.isna(
        )]
        data_dfs[0] = data_dfs[0][~data_dfs[0].Peptide.isin(
            low_complex_targets)].reset_index(drop=True)

        #doing peptide level competition
        df_all = pf.peptide_level(
            data_dfs[0].copy(), peptide_list_dfs[0].copy(), remove)

        PSMs = data_dfs[0].copy()
        PSMs['rank'] = PSMs['SpecId'].apply(
            lambda x: int(x[-1]))
        PSMs = PSMs[PSMs['rank'] == 1].reset_index(drop=True)
        PSMs.drop(remove, axis=1, inplace=True, errors='ignore')

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

        #do SVM
        train_power, std_power, true_power, df_new, train_all_new, model, columns_trained = pf.do_svm(df_all_scale.copy(), train_all.copy(), df_all.copy(), folds=folds, Cs=[
            0.1, 1, 10], p=p_init, total_iter=total_iter, alpha=FDR_threshold, train_alpha=train_FDR_threshold)

        df_new = df_new.loc[df_new.q_val <= FDR_threshold]

    else:
        sys.stderr.write("Reading in search file(s). \n")
        logging.info("Reading in search file(s).")
        data_dfs = []
        for search_file in search_files:
            data_df = uf.read_pin(search_file)
            #removing flanking aa
            data_df['Peptide'] = data_df['Peptide'].apply(
                lambda x: x[2:(len(x) - 2)])
            data_dfs.append(data_df)

        #averaging the two tailor PSMs.
        data_dfs[0] = data_dfs[0].merge(
            data_dfs[1][['SpecId', 'TailorScore']], how='left', on='SpecId')
        
        sys.stderr.write("Fixing the Tailor Scores. \n")
        logging.info("Fixing the Tailor Scores.")
        
        data_dfs[0]['TailorScore'] = (
            data_dfs[0].TailorScore_x + data_dfs[0].TailorScore_y)/2
        data_dfs[0].drop(['TailorScore_x', 'TailorScore_y'],
                         inplace=True, axis=1)
        data_dfs.pop(1)
        
        #do PSM level competition involving all PSMs
        data_df = pf.PSM_level(data_dfs[0].copy(), data_dfs[1:].copy(), top=1)
        PSMs = data_df.copy()
        
        sys.stderr.write("Reading in peptide list(s). \n")
        logging.info("Reading in peptide list(s).")
        peptide_list_dfs = []
        for td_list in td_lists:
            peptide_list_df = pd.read_table(td_list)
            if 'decoy(s)' in peptide_list_df.columns:
                peptide_list_df.rename(columns={'decoy(s)': 'decoy'}, inplace=True)
            peptide_list_df.drop_duplicates(inplace=True)
            peptide_list_dfs.append(peptide_list_df)
        
        #removing low complexity targets that do not produce a decoy
        low_complex_targets = peptide_list_dfs[0].target[peptide_list_dfs[0].decoy.isna(
        )]
        data_dfs[0] = data_dfs[0][~data_dfs[0].Peptide.isin(
            low_complex_targets)].reset_index(drop=True)

        

        PSMs = data_dfs[0].copy()
        PSMs['rank'] = PSMs['SpecId'].apply(
            lambda x: int(x[-1]))
        PSMs = PSMs[PSMs['rank'] == 1].reset_index(drop=True)
        PSMs.drop(remove, axis=1, inplace=True, errors='ignore')

        #doing multi peptide level competition
        df_all = pf.peptide_level(
            data_df.copy(), peptide_list_dfs.copy(), remove)

        #do scale
        df_all_scale, scale = pf.do_scale(df_all.copy())

        if p_init_default:
            p_init = 0.5

        #create decoys at pseudolevel
        rand_indxs = np.random.choice([True, False], replace = True, size = sum((df_all_scale['Label']
                                      == -1)), p = [p_init, 1 - p_init])
        
        train_all = df_all_scale.loc[(df_all_scale['Label']
                                      == -1)].copy()
        train_all = train_all.loc[rand_indxs].copy()

        #do svm
        train_power, std_power, true_power, df_new, train_all_new, model, columns_trained = pf.do_svm(df_all_scale.copy(), train_all.copy(), df_all.copy(), folds=folds, Cs=[
            0.1, 1, 10], p=p_init, total_iter=total_iter, alpha=FDR_threshold, train_alpha=train_FDR_threshold, mult=2)

        df_new = df_new.loc[df_new.q_val <= FDR_threshold]

    if df_new.shape[0] > 0 and get_psms:
        sys.stderr.write(
            "Reporting all PSMs within each mass-cluster associated to a discovered peptide. \n")
        logging.info(
            "Reporting all PSMs within each mass-cluster associated to a discovered peptide. ")
        originally_discovered = df_new.loc[df_new.Label == 1, 'Peptide'].copy()
        #originally_discovered = originally_discovered.str.replace(
        #    "\\[|\\]|\\.|\\d+", "", regex=True)
        df_extra = uf.create_cluster(PSMs.copy(
        ), scale, originally_discovered, model, isolation_window, columns_trained)
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

        #write results
        if output_dir != '.':
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
            df_new[df_new.Label == 1].to_csv(output_dir + "/" + file_root +
                                             ".peptides.txt", header=True, index=False, sep='\t')
            df_final[df_final.Label == 1].to_csv(output_dir + "/" + file_root +
                                                 ".psms.txt", header=True, index=False, sep='\t')

        else:
            df_new.to_csv(output_dir + "/" + file_root +
                          ".peptides.txt", header=True, index=False, sep='\t')
            df_final[df_final.Label == 1].to_csv(output_dir + "/" + file_root +
                                                 ".psms.txt", header=True, index=False, sep='\t')
    else:

        #write results
        if output_dir != '.':
            if not os.path.isdir(output_dir):
                os.mkdir(output_dir)
            df_new[df_new.Label == 1].to_csv(output_dir + "/" + file_root +
                                             ".peptides.txt", header=True, index=False, sep='\t')

        else:
            df_new.to_csv(output_dir + "/" + file_root +
                          ".peptides.txt", header=True, index=False, sep='\t')

    end_time = time.time()

    logging.info("Elapsed time: " +
                 str(round(end_time - start_time, 2)) + " s")
    sys.stderr.write("Elapsed time: " +
                     str(round(end_time - start_time, 2)) + " s \n")


if __name__ == "__main__":
    main()

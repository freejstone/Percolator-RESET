#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 18:27:26 2023

@author: jackfreestone

This module performs the percolator algorithm with strict FDR control using just the open-search for now
"""

import os
import time
import sys
import numpy as np
import pandas as pd
import random
import logging
import utility_functions as uf
import super_percolator_open_functions as spf

USAGE = """USAGE: python3 do_super_percolator.py [options] <wide> <matching>

    To be added.
            
"""

#########################################################################################################


def main():
    global USAGE

    start_time = time.time()

    # Set default values for parameters.
    FDR_threshold = 0.01
    folds = 3
    total_iter = 10
    kernel = 'linear'
    train_FDR_threshold = 0.01
    degree = None
    account_mods = True
    precursor_bin_width = 1.0005079/4
    neighbour_remove = True
    remove = ['bins', 'pi_0', 'freq', 'min_tailor_score', 'min_xcorr_score']
    thresh = 0.05
    top_positive = True
    n_processes = 1
    output_dir = '.'
    file_root = 'super_percolator'
    static_mods = {'C': 57.02146}
    overwrite = False
    seed = None
    command_line = ' '.join(sys.argv)
    stratified = False
    iterate = False
    FDR_grid = [0.25, 0.1]
    p_init = 0.25

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
        elif (next_arg == "--top_positive"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                top_positive = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                top_positive = False
            else:
                sys.stderr.write("Invalid argument for --top_positive")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--kernel"):
            kernel = str(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--train_FDR_threshold"):
            train_FDR_threshold = float(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--degree"):
            degree = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--account_mods"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                account_mods = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                account_mods = False
            else:
                sys.stderr.write("Invalid argument for --account_mods")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--precursor_bin_width"):
            precursor_bin_width = float(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--neighbour_remove"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                neighbour_remove = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                neighbour_remove = False
            else:
                sys.stderr.write("Invalid argument for --neighbour_remove")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--remove"):
            remove = str(sys.argv[0]).split(',')
            sys.argv = sys.argv[1:]
        elif (next_arg == "--n_processes"):
            n_processes = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--output_dir"):
            output_dir = str(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--file_root"):
            file_root = str(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == "--static_mods"):
            static_mods = uf.parse_static_mods(sys.argv[0])
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
        elif (next_arg == "--stratified"):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                stratified = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                stratified = False
            else:
                sys.stderr.write("Invalid argument for --stratified")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == '--seed'):
            seed = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == '--stratified'):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                stratified = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                stratified = False
            else:
                sys.stderr.write("Invalid argument for --stratified")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == '--iterate'):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                iterate = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                iterate = False
            else:
                sys.stderr.write("Invalid argument for --iterate")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == '--FDR_grid'):
            FDR_grid = [float(j) for j in sys.argv[0].split(',')]
            sys.argv = sys.argv[1:]
        elif (next_arg == "--p_init"):
            p_init = float(sys.argv[0])
            sys.argv = sys.argv[1:]
        else:
            sys.stderr.write("Invalid option (%s)" % next_arg)
            sys.exit(1)
    if (len(sys.argv) == 2):
        search_file_open = sys.argv[0]
        td_list = sys.argv[1]
    else:
        #sys.stderr.write('Version: ' + str(__version__) + " \n")
        sys.stderr.write(USAGE)
        sys.exit(1)

    #setting seed for reproducibility
    if type(seed) == int:
        random.seed(seed)
        np.random.seed(seed)

    #print meta information, checking directory and printing warnings
    uf.print_info(command_line, output_dir, file_root, overwrite,
                  account_mods, None, search_file_open)

    sys.stderr.write("Reading in search file and peptide list. \n")
    logging.info("Reading in search file and peptide list.")

    open_df = uf.read_pin(search_file_open)  # reading
    peptide_list_df = pd.read_table(td_list)  # reading

    if 'decoy(s)' in peptide_list_df.columns:
        peptide_list_df.rename(columns={'decoy(s)': 'decoy'}, inplace=True)

    peptide_list_df.drop_duplicates(inplace=True)

    #removing flanking aa
    open_df['Peptide'] = open_df['Peptide'].apply(
        lambda x: x[2:(len(x) - 2)])

    #doing peptide level competition and get bins/freq
    df_all = spf.peptide_level(
        open_df.copy(), peptide_list_df.copy(), precursor_bin_width=precursor_bin_width)
    
    if stratified:
        #get second sample proportion
        p_next = 1 - 1/(2*(1 - p_init))

        #do scale
        df_all = spf.do_scale(df_all)

        q_val = uf.TDC_flex_c(
            df_all.Label == -1, df_all.Label == 1)

        df_below = df_all[~(q_val <= FDR_grid[0])].copy()

        df_above = df_all[(q_val <= FDR_grid[0])].copy()

        #create target-decoys at pseudolevel
        train_decoys_above = df_above.loc[(df_above['Label']
                                           == -1)].sample(frac=p_init).copy()

        train_below = df_below

        train_all = pd.concat([train_decoys_above, train_below])

        train_power, std_power, true_power, df_new, train_all_new = spf.do_svm(df_above.copy(), train_all.copy(), folds=folds, Cs=[
            0.1, 1, 10], p=p_init, total_iter=total_iter, kernel=kernel, alpha=FDR_threshold, train_alpha=train_FDR_threshold, degree=degree, remove=remove, top_positive=False)
        #df_new, train_decoys_new = do_lda(df_above, train_all, total_iter=3, p=0.5, alpha=0.01, train_alpha=0.01, remove=remove, top_positive=False, qda=False)
        
        df_new['trained'] = 0  # has not been trained
        train_all_new['trained'] = 1  # has been trained
        df_new = pd.concat([df_new, train_all_new])
        df_new = spf.do_scale(df_new)

        df_new = df_new.sort_values(
            by='SVM_score', ascending=False).reset_index(drop=True)

        train_all_new = df_new[df_new.trained == 1].copy()
        df_new = df_new[df_new.trained == 0].copy()

        q_val = uf.TDC_flex_c(
            df_new.Label == -1, df_new.Label == 1, c=1/(2 - p_init), lam=1/(2 - p_init))

        train_all = pd.concat([train_all_new, df_new[q_val > FDR_grid[1]].copy(
        ), df_new[(q_val <= FDR_grid[1]) & (df_new.Label == -1)].sample(frac=p_next).copy()])

        train_power_next, std_power_next, true_power_next, df_new, train_all_new = spf.do_svm(df_new.copy(), train_all.copy(), folds=folds, Cs=[
            0.1, 1, 10], p=p_next, total_iter=total_iter, kernel=kernel, alpha=FDR_threshold, train_alpha=train_FDR_threshold, degree=degree, remove=remove, top_positive=False)
        
        df_new['trained'] = 0  # has not been trained
        train_all_new['trained'] = 1  # has been trained
        df_new = pd.concat([df_new, train_all_new])
        df_new = spf.do_scale(df_new)

        df_new = df_new.sort_values(
            by='SVM_score', ascending=False).reset_index(drop=True)
        
        if iterate:
            R = sum((df_new.Label == 1) & (df_new.trained == 0))
            V = sum((df_new.Label == -1) & (df_new.trained == 0))
            c = 1/(2 - p_next)
            lam = 1/(2 - p_next)
            est_FDP = (V + 1)*c/((1 - lam)*max(R, 1))
            
            while est_FDP > 0.01:
                indxs = df_new[(df_new['trained'] == 0)].index[-2:]
    
                spec_ids = df_new.SpecId.loc[indxs]
    
                train_all = pd.concat([df_new[df_new.trained == 1].copy(), df_new[df_new.SpecId.isin(spec_ids)]])
    
                df_new = df_new.drop('trained', axis = 1)
                train_all = train_all.drop('trained', axis = 1)
    
                train_power_next, std_power_next, true_power_next, df_new, train_all = spf.do_svm(df_new.copy(), train_all.copy(), folds=folds, Cs=[
                            0.1, 1, 10], p = p_next, total_iter=10, kernel=kernel, alpha=FDR_threshold, train_alpha=train_FDR_threshold, degree=degree, remove=remove, top_positive=False)
    
                #df_new, train_decoys_new = do_lda(df_new_above, train_decoys_new, total_iter=3, p=0.5, alpha=0.01, train_alpha=0.01, remove=remove, top_positive=False, qda=False)
    
                df_new['trained'] = 0 #has not been trained
                train_all_new['trained'] = 1 #has been trained
                df_new = pd.concat([df_new, train_all_new])
                df_new = spf.do_scale(df_new)
    
                df_new = df_new.sort_values(
                    by='SVM_score', ascending=False).reset_index(drop=True)
    
                R = sum((df_new.Label == 1) & (df_new.trained == 0))
                V = sum((df_new.Label == -1) & (df_new.trained == 0))
                c = 1/(2 - p_next)
                lam = 1/(2 - p_next)
                est_FDP = (V + 1)*c/((1 - lam)*max(R, 1))
                
                sys.stderr.write("estimated FDP is %s. \n" %(est_FDP))
                logging.info("estimated FDP is %s." %(est_FDP))
            
            p_init = p_next
        
        df_new = df_new[df_new.trained == 0].copy()
        
        q_val = uf.TDC_flex_c(
            df_new.Label == -1, df_new.Label == 1, c=1/(2 - p_init), lam=1/(2 - p_init))

        df_new['q_val'] = q_val
        
    else:
        
        #do scale
        df_all = spf.do_scale(df_all)

        q_val = uf.TDC_flex_c(
            df_all.Label == -1, df_all.Label == 1)

        #create target-decoys at pseudolevel
        train_all = df_all.loc[(df_all['Label']
                                           == -1)].sample(frac=p_init).copy()

        train_power, std_power, true_power, df_new, train_all_new = spf.do_svm(df_all.copy(), train_all.copy(), folds=folds, Cs=[
            0.1, 1, 10], p=p_init, total_iter=10, kernel=kernel, alpha=FDR_threshold, train_alpha=train_FDR_threshold, degree=degree, remove=remove, top_positive=False)
        
        df_new = df_new.sort_values(
            by='SVM_score', ascending=False).reset_index(drop=True)
        
        q_val = uf.TDC_flex_c(
            df_new.Label == -1, df_new.Label == 1, c=1/(2 - p_init), lam=1/(2 - p_init))
        
        df_new['q_val'] = q_val
    

    obs_power = sum((q_val <= FDR_threshold) & (df_new.Label == 1))

    sys.stderr.write("%s peptides discovered at %s FDR. \n" %(obs_power, FDR_threshold))
    logging.info("%s peptides discovered at %s FDR." %(obs_power, FDR_threshold))

    #write results
    if output_dir != './':
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        df_new.to_csv(output_dir + "/" + file_root +
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
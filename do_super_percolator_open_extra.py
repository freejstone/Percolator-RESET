#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 12:34:20 2023

@author: jfre0619

This module performs the percolator algorithm with strict FDR control using just the open-search for now with extra decoy
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
    psm_level = True
    p_init = 0.5
    open_narrow = 'open'
    keep_hidden = True

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
        elif (next_arg == '--p_init'):
            p_init = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == '--psm_level'):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                psm_level = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                psm_level = False
            else:
                sys.stderr.write("Invalid argument for --psm_level")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        elif (next_arg == "--open_narrow"):
            if str(sys.argv[0]) in ['open', 'narrow']:
                open_narrow = str(sys.argv[0])
            else:
                sys.stderr.write("Invalid argument for --open_narrow")
            sys.argv = sys.argv[1:]
        elif (next_arg == '--keep_hidden'):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                keep_hidden = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                keep_hidden = False
            else:
                sys.stderr.write("Invalid argument for --keep_hidden")
                sys.exit(1)
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

    #print meta information, checking directory and printing warnings
    if open_narrow == 'open':
        uf.print_info(command_line, output_dir, file_root, overwrite,
                      account_mods, None, search_files)
    elif open_narrow == 'narrow':
        uf.print_info(command_line, output_dir, file_root, overwrite,
                      account_mods, search_files, None)


    sys.stderr.write("Reading in search file and peptide list. \n")
    logging.info("Reading in search file and peptide list.")
    
    data_dfs = []
    peptide_list_dfs = []
    for search_file in search_files:
        data_dfs.append(uf.read_pin(search_file))
    for td_list in td_lists:
        peptide_list_df = pd.read_table(td_list)
        if 'decoy(s)' in peptide_list_df.columns:
            peptide_list_df.rename(columns = {'decoy(s)':'decoy'}, inplace = True)
        peptide_list_df.drop_duplicates(inplace = True)
        peptide_list_dfs.append(peptide_list_df)
    
    if psm_level:
        
        #do PSM level competition
        if keep_hidden:
            data_df = spf.PSM_level(data_dfs[0], data_dfs[1], top = 1)
            
            #removing flanking aa
            data_df['Peptide'] = data_df['Peptide'].apply(
                lambda x: x[2:(len(x) - 2)])
            
            df_extra_decoy = data_dfs[2].copy()
            
            #remove repeated decoy peptides
            df_extra_decoy['Peptide'] = df_extra_decoy['Peptide'].apply(
                lambda x: x[2:(len(x) - 2)])
            df_extra_decoy = df_extra_decoy[~(df_extra_decoy.Peptide.isin(peptide_list_dfs[0].decoy))]
            df_extra_decoy = spf.PSM_level(df_extra_decoy.copy(), None, top = 1) #this just gets the top 1 PSM
            
            #get pseudo PSM level competition
            df_pseudo, df_all = spf.pseudo_PSM_level(data_df.copy(), df_extra_decoy.copy(), 1, precursor_bin_width)
            
        if not keep_hidden:
            data_df = spf.PSM_level(data_dfs[0], data_dfs[1:], top = 1)
            
            #removing flanking aa
            data_df['Peptide'] = data_df['Peptide'].apply(
                lambda x: x[2:(len(x) - 2)])
            
            df_all = data_df[data_df.filename == 0].copy().reset_index(drop = True)
            
            df_pseudo = data_df.copy()
            df_pseudo.loc[(df_pseudo.filename == 0),'Label'] = 1
            df_pseudo.loc[(df_pseudo.filename == 1),'Label'] = -1
        
        df_pseudo, df_all = spf.do_scale(df_pseudo.copy(), df_all.copy())
        
        #Doing iterative PSM
        train_power, std_power, true_power, real_df = spf.do_svm_extra(df_pseudo.copy(), df_all.copy(), folds=folds, Cs=[
            0.1, 1, 10], total_iter=total_iter, kernel=kernel, alpha=FDR_threshold, train_alpha=train_FDR_threshold, degree=degree, remove=remove, top_positive=False)
        
        #doing peptide-level competition now
        real_df_peptide = spf.peptide_level(real_df.copy(), peptide_list_dfs[0].copy(), precursor_bin_width=precursor_bin_width, open_narrow = open_narrow)
        
        #the new direction
        real_df_peptide.sort_values(
            by='SVM_score', ascending=False).reset_index(drop=True)
        
        #get_qvals
        q_val = uf.TDC_flex_c(
            real_df_peptide.Label == -1, real_df_peptide.Label == 1, c=1/2, lam=1/2)
        
        real_df_peptide['q_vals'] = q_val
        
        obs_power = sum((q_val <= FDR_threshold) & (real_df_peptide.Label == 1))

        sys.stderr.write("%s peptides discovered at %s FDR. \n" %(obs_power, FDR_threshold))
        logging.info("%s peptides discovered at %s FDR." %(obs_power, FDR_threshold))

    else:
        
        #do PSM level competition involving all PSMs
        data_dfs = spf.PSM_level(data_dfs[0], data_dfs[1:], top = 1)
        
        #removing flanking aa
        data_dfs['Peptide'] = data_dfs['Peptide'].apply(
            lambda x: x[2:(len(x) - 2)])
        
        #doing multi peptide level competition and get bins/freq for al
        df_all = spf.peptide_level(
            data_dfs.copy(), peptide_list_dfs.copy(), precursor_bin_width=precursor_bin_width, open_narrow = open_narrow)
        
        #do scale
        df_all = spf.do_scale(df_all)

        #create target-decoys at pseudolevel
        train_all = df_all.loc[(df_all['Label']
                                           == -1)].sample(frac=p_init).copy()
        
        #do svm with
        train_power, std_power, true_power, real_df_peptide, train_all_new = spf.do_svm(df_all.copy(), train_all.copy(), folds=folds, Cs=[
            0.1, 1, 10], p=p_init, total_iter=10, kernel=kernel, alpha=FDR_threshold, train_alpha=train_FDR_threshold, degree=degree, remove=remove, top_positive=False, mult = 2)
        
        real_df_peptide = real_df_peptide.sort_values(
            by='SVM_score', ascending=False).reset_index(drop=True)
        
        q_val = uf.TDC_flex_c(
            real_df_peptide.Label == -1, real_df_peptide.Label == 1, c=1/(2*(1 - p_init) + 1), lam=1/(2*(1 - p_init) + 1))
        
        real_df_peptide['q_val'] = q_val
        
   
    #write results
    if output_dir != './':
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        real_df_peptide.to_csv(output_dir + "/" + file_root +
                           ".peptides.txt", header=True, index=False, sep='\t')

    else:
        real_df_peptide.to_csv(output_dir + "/" + file_root +
                           ".peptides.txt", header=True, index=False, sep='\t')

    end_time = time.time()

    logging.info("Elapsed time: " +
                 str(round(end_time - start_time, 2)) + " s")
    sys.stderr.write("Elapsed time: " +
                     str(round(end_time - start_time, 2)) + " s \n")


if __name__ == "__main__":
    main()
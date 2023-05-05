#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 12:30:14 2023

@author: jfre0619
"""
import os
import time
import sys
import numpy as np
import pandas as pd
import random
import logging
import utility_functions as uf
import super_percolator_functions_extra as spf

USAGE = """USAGE: python3 do_super_percolator.py [options] <narrow> <wide> <matching>

  This script implements the super percolator algorithm. The first input file
  is the narrow search file output, the second input file is the open search file
  output, and the last input file contains the target-decoy peptide
  pairs. Output is a list of peptides discovered at a user-specified FDR-level.

  Options:
      
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
    remove = ['bins', 'pi_0', 'enzInt']
    thresh = 0.05
    top_positive = True
    n_processes = 1
    output_dir = '.'
    file_root = 'super_percolator'
    static_mods = {'C': 57.02146}
    overwrite = False
    seed = None
    command_line = ' '.join(sys.argv)
    svm = True

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
        elif (next_arg == '--seed'):
            seed = int(sys.argv[0])
            sys.argv = sys.argv[1:]
        elif (next_arg == '--svm'):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                svm = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                svm = False
            else:
                sys.stderr.write("Invalid argument for --svm")
                sys.exit(1)
            sys.argv = sys.argv[1:]
        else:
            sys.stderr.write("Invalid option (%s)" % next_arg)
            sys.exit(1)
    if (len(sys.argv) == 3):
        search_files_narrow = str(sys.argv[0]).split(',')
        search_files_open = str(sys.argv[1]).split(',')
        td_lists = str(sys.argv[2]).split(',')
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
                  account_mods, search_files_narrow, search_files_open)

    sys.stderr.write("Reading in search files and peptide list. \n")
    logging.info("Reading in search files and peptide list.")
    
    #reading
    narrow_dfs = []
    open_dfs = []
    peptide_list_dfs = []
    for search_file_narrow in search_files_narrow:
        narrow_dfs.append(uf.read_pin(search_file_narrow))
    for search_file_open in search_files_open:
        open_dfs.append(uf.read_pin(search_file_open))
    for td_list in td_lists:
        peptide_list_df = pd.read_table(td_list)
        if 'decoy(s)' in peptide_list_df.columns:
            peptide_list_df.rename(columns = {'decoy(s)':'decoy'}, inplace = True)
        peptide_list_df.drop_duplicates(inplace = True)
        peptide_list_dfs.append(peptide_list_df)
    
    #do PSM_level competitions
    narrow_df_comp = spf.PSM_level(narrow_dfs[0], narrow_dfs[1], top = 1)
    open_df_comp = spf.PSM_level(open_dfs[0], open_dfs[1], top = 5)
    
    #doing filtering
    df_all = uf.filter_narrow_open(narrow_df_comp, open_df_comp, thresh, n_processes,
                                   neighbour_remove, tide_used='tide', static_mods=static_mods)
    
    #doing filtering using extra decoy
    df_extra_decoy_narrow, df_extra_decoy_open = narrow_dfs[-1].copy(), open_dfs[-1].copy()
    
    df_extra_decoy_narrow = spf.get_rank(df_extra_decoy_narrow, 1) #getting ranks
    df_extra_decoy_open = spf.get_rank(df_extra_decoy_open, 1) 
    
    df_extra_decoy_narrow['n_o'] = 1 #getting narrow_open label
    df_extra_decoy_open['n_o'] = 0
    
    df_extra_decoy = pd.concat([df_extra_decoy_narrow, df_extra_decoy_open]) #combining narrow_open extra decoys
    
    #remove repeated decoy peptides
    df_extra_decoy['Peptide'] = df_extra_decoy['Peptide'].apply(
        lambda x: x[2:(len(x) - 2)])
    df_extra_decoy = df_extra_decoy[~(df_extra_decoy.Peptide.isin(peptide_list_dfs[0].decoy))]
    
    df_extra_decoy.reset_index(drop = True, inplace = True)
        
    #get pseudo PSM level competition
    df_pseudo, df_all = spf.pseudo_PSM_level(df_all.copy(), df_extra_decoy.copy(), 2, precursor_bin_width)
                                 
    #Doing iterative PSM
    train_power, std_power, true_power, real_df = spf.do_iterative_svm_cv(df_pseudo.copy(), df_all.copy(), folds=folds, Cs=[
        0.1, 1, 10], total_iter=total_iter, kernel=kernel, alpha=FDR_threshold, train_alpha=train_FDR_threshold, degree=degree, remove=remove, top_positive=top_positive)
    
    #doing peptide-level competition now
    real_df_peptide = spf.peptide_level(real_df.copy(), peptide_list_dfs[0].copy())
    
    #the new direction
    real_df_peptide.sort_values(
        by='SVM_scores', ascending=False).reset_index(drop=True)
    
    #get_qvals
    q_val = uf.TDC_flex_c(
        real_df_peptide.Label == -1, real_df_peptide.Label == 1, c=1/2, lam=1/2)
    
    real_df_peptide['q_vals'] = q_val
    
    obs_power = sum((q_val <= FDR_threshold) & (real_df_peptide.Label == 1))

    sys.stderr.write("%s peptides discovered at %s FDR. \n" %(obs_power, FDR_threshold))
    logging.info("%s peptides discovered at %s FDR." %(obs_power, FDR_threshold))

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

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
import re
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
    remove = ['bins', 'pi_0', 'min_tailor_score', 'min_xcorr_score']
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
    df_extra_decoy = uf.filter_narrow_open(narrow_dfs[-1], open_dfs[-1], thresh, n_processes,
                                   neighbour_remove, tide_used='tide', static_mods=static_mods)
    
    #get extra features here like bin, freq, pi_0
    df_pseudo_target = df_all.copy()
    df_pseudo_target = df_pseudo_target[df['rank'] == 1]
    df_pseudo_target['Label'] = 1
    df_extra_decoy = 
    
    
    #do pseudo-competition at the PSM level
    
    #do iterative SVM
    
    #do peptide-level competition
    
    #do TDC
    
    #need to remove df_extra_decoy that overlap with one of them
    #df_extra_decoy = df_extra_decoy[~(df_extra_decoy.Peptide.isin(peptide_list_dfs[0].decoy))].reset_index(drop = True)

    #doing peptide level competition base
    df_all = peptide_level(df_all, peptide_list_dfs[0]).reset_index(drop = True)
    
    #doing the same for the extra decoy -- this is really just to remove low ranked peptides
    df_extra_decoy = peptide_level(df_extra_decoy, peptide_list_dfs[-1]).reset_index(drop = True)
    
    #preparing pseudo level competition
    df_all_pseudo = df_all.copy().drop(['min_tailor_score', 'min_xcorr_score'], axis = 1)
    df_all_pseudo['Label'] = 1
    df_all_pseudo = df_all_pseudo[df_all_pseudo.original_target.isin(df_extra_decoy.original_target)]
    df_extra_decoy = df_extra_decoy[df_extra_decoy.original_target.isin(df_all_pseudo.original_target)]
    df_all_pseudo = pd.concat([df_all_pseudo, df_extra_decoy])
    
    #getting bin, freq, pi0
    df_all_pseudo = spf.get_bin_freq_pi0(df_all_pseudo, precursor_bin_width = precursor_bin_width)
    df_all[['bins', 'freq', 'pi_0']] = df_all_pseudo[['bins', 'freq', 'pi_0']][df_all_pseudo.Label == 1]
    
    #doing pseudo level competition
    df_all_pseudo = peptide_level(df_all_pseudo, peptide_list_dfs[-1])
    
    #drop original_target
    df_all_pseudo.drop('original_target', axis = 1, inplace = True)
    df_all.drop('original_target', axis = 1, inplace = True)
    
    df_all_pseudo.drop('n_o', axis = 1, inplace = True)
    df_all.drop('n_o', axis = 1, inplace = True)
    
    
    train_power, std_power, true_power, discoveries = do_iterative_svm_cv(df_all_pseudo.copy(), df_all.copy(), folds=folds, Cs=[
        0.1, 1, 10], total_iter=total_iter, kernel=kernel, alpha=FDR_threshold, train_alpha=train_FDR_threshold, degree=degree, remove=remove, top_positive=top_positive)
    power = pd.DataFrame(zip(train_power, std_power, true_power), columns=[
                         'train_power', 'std_power', 'true_power'])
        

    #write results
    if output_dir != './':
        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)
        discoveries.to_csv(output_dir + "/" + file_root +
                           ".peptides.txt", header=True, index=False, sep='\t')
        power.to_csv(output_dir + "/" + file_root + ".power.txt",
                     header=True, index=False, sep='\t')

    else:
        discoveries.to_csv(output_dir + "/" + file_root +
                           ".peptides.txt", header=True, index=False, sep='\t')

    end_time = time.time()

    logging.info("Elapsed time: " +
                 str(round(end_time - start_time, 2)) + " s")
    sys.stderr.write("Elapsed time: " +
                     str(round(end_time - start_time, 2)) + " s \n")


if __name__ == "__main__":
    main()

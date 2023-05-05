#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 17:15:01 2023

@author: jackfreestone

This module performs the percolator algorithm with strict FDR control 
"""
import os
import time
import sys
import numpy as np
import pandas as pd
import random
import logging
import utility_functions as uf
import super_percolator_functions as spf

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
    stratified = False

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
        elif (next_arg == '--svm'):
            if str(sys.argv[0]) in ['t', 'T', 'true', 'True']:
                svm = True
            elif str(sys.argv[0]) in ['f', 'F', 'false', 'False']:
                svm = False
            else:
                sys.stderr.write("Invalid argument for --svm")
                sys.exit(1)
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
        else:
            sys.stderr.write("Invalid option (%s)" % next_arg)
            sys.exit(1)
    if (len(sys.argv) == 3):
        search_file_narrow = sys.argv[0]
        search_file_open = sys.argv[1]
        td_list = sys.argv[2]
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
                  account_mods, search_file_narrow, search_file_open)

    sys.stderr.write("Reading in search files and peptide list. \n")
    logging.info("Reading in search files and peptide list.")

    narrow_df = uf.read_pin(search_file_narrow)  # reading
    open_df = uf.read_pin(search_file_open)  # reading
    peptide_list_df = pd.read_table(td_list)  # reading
    
    if 'decoy(s)' in peptide_list_df.columns:
        peptide_list_df.rename(columns = {'decoy(s)':'decoy'}, inplace = True)
    
    peptide_list_df.drop_duplicates(inplace = True)
    
    #doing filtering
    df_all = uf.filter_narrow_open(narrow_df, open_df, thresh, n_processes,
                                   neighbour_remove, tide_used='tide', static_mods=static_mods)

    #doing peptide level competition
    df_all = spf.peptide_level(
        df_all, peptide_list_df, precursor_bin_width=precursor_bin_width)
    
    if svm:
        if stratified:
            train_power, std_power, true_power, discoveries = spf.do_iterative_svm_cv_stratified(df_all, folds=folds, Cs=[
                0.1, 1, 10], p = 0.5, total_iter=total_iter, kernel=kernel, alpha=FDR_threshold, train_alpha=train_FDR_threshold, degree=degree, remove=remove, top_positive=top_positive)
        else:
            train_power, std_power, true_power, discoveries = spf.do_iterative_svm_cv(df_all, folds=folds, Cs=[
                0.1, 1, 10], p = 0.5, total_iter=total_iter, kernel=kernel, alpha=FDR_threshold, train_alpha=train_FDR_threshold, degree=degree, remove=remove, top_positive=top_positive)
        power = pd.DataFrame(zip(train_power, std_power, true_power), columns=[
                             'train_power', 'std_power', 'true_power'])
    else:
        true_power, discoveries = spf.do_iterative_lda_cv(df_all, total_iter=total_iter, p = 0.5, alpha=FDR_threshold, train_alpha=train_FDR_threshold, remove=remove, top_positive=top_positive, qda = True)
        power = pd.DataFrame(true_power, columns=[
                             'true_power'])

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

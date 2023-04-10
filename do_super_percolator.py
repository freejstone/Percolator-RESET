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
      
    --output_dir <string> The file-path of the output
                          Default = './'.
                          
    --file_root <string>  The file prefix of the output files
                          Default = 'super_percolator'.

    --FDR_threshold <value>     A user-specified threshold for which the reported
                                peptides will have FDR below this threshold.
                                Default = 0.01.

    --K <integer>         The number of recently observed peptides used
                          to estimate the probability that the next
                          peptide is a target or decoy.
                          Default = 40.

    --tops_gw <integer>   The number of top PSMs for each scan in the open 
                          search that will be used by group-walk.
                          Default = 2.
    
    --tops_open <integer> The number of top PSMs in the open search used in
                          the neighbour-filtering process.
                          Default = 5.                         
                          
    --score <string>      Either 'tailor_score', 'xcorr_score', 'e-value' or
                          'hyperscore'. The score that will be used in the 
                          peptide-level competition and subsequent Group 
                          construction and Group-walk algorithm. If 
                          'tailor_score', it is assumed the search files are 
                          derived from Tide. If 'xcorr_score', either Tide 
                          search or Comet is assumed to be used. If 'e-value',
                          Comet is assumed. If 'hyperscore' it is assumed the
                          search files are derived from MS-Fragger.
                          Default = 'tailor_score'.
    
    --account_mods <T|F>  To determine whether the group-walk algorithm
                          selects the best PSM among the equivalent
                          classes of peptides which are equal up to
                          variable modification, or not.
                          Default = T.
                          
    --isolation_window <value>    The left and right isolation window offsets
                                    used in tandem MS/MS. Values should be comma-
                                    separated.
                                    Default = 2,2
                                    
                                    precursors according to their mass
                                    for subsequent competition. The first value should
                                    be the left-offset of the isolation window and the
                                    second value should be the right-offset. The
                                    two values should be comma-separated.
                                    Default = 2,2.
                          
    --precursor_bin_width <value>   To determine the size of the bins
                                    used to discretize the mass-
                                    differences between the sample
                                    and theoretical spectra.
                                    Default = 1.0005079/4.
                          
    --print_chimera <T|F>          To determine whether we print the number
                                   of scans that have more than 1 peptide
                                   discovered at the 1% and 5% FDR level
                                   to the log file.
                                   Default = T.
            
    --group_thresh <value>         The p-value threshold used to determine
                                   whether groups are merged or kept
                                   separate in the KS test.
                                   Default = 0.01.
                                   
    --print_group_pi0 <T|F>        To determine whether the group-level
                                   proportion of pi_0 is printed to the 
                                   log file.
                                   Default = T.
    
    --min_group_size <integer>     The number of multiples of K that is
                                   used to determine the minimum size of
                                   each group. See option --K.
                                   Default = 2.
    
    --neighbour_remove <T|F>       If true, for each scan, we successively
                                   move down the list of PSMs associated with
                                   each scan, ordered in terms of the score
                                   from highest to lowest, and compute
                                   a similarity score between the current
                                   peptide and the previous peptide(s).
                                   If one of the similarity score(s) exceeds
                                   a certain threshold, the PSM associated 
                                   with the current peptide is thrown out,
                                   and we proceed to the next.
                                   Default = T.
    
    --thresh <value>      The similarity score used as a threshold to filter
                          out neighbouring peptides.
                          Default = 0.05.
                          
    --return_filt_search <T|F>  Whether or not to return filtered narrow
                                and open search files.
                                Default = F.
                                
    --return_frontier <T|F>     The sequence of indices describing the
                                positions of the frontier used by Groupwalk
                                is returned as a .txt file to the output
                                directory.
                                Default = F.
    
    --n_processes <integer>     The number of threads, used in the filtering
                                process.
                                Default = 1.
                                
    --static_mods <string>      Of the form X:[+-]A where X is the amino acid,
                                or rather "cterm" or "nterm" if it is a
                                modification on the C-terminal or N-terminal.
                                A is the absolute mass shift in Daltons.
                                [+-] indicates whether the mass shift is
                                positive or negative. C+57.02146 is always
                                included by default. Variable modifications
                                do not need specification (they are accounted
                                for via the search files). List mods in comma
                                separated format, e.g.
                                nterm:10,cterm:-20,L:50.
                                Default = None.
                                
    --return_mass_mod_hist <T|F>   A matplotlib histogram of mass 
                                   modifications is returned to the user's
                                   directory.
                                   Default = F.
                                
    --dcy_prefix <string>       The prefix used for the decoy proteins.
                                Default = 'decoy_'

    --return_decoys <T|F>       Also report decoys used to estimate the 
                                number of false discoveries. 
                                Default = F. 
                                
    --overwrite <T|F>     Gives option to overwrite existing files in
                          directory.
                          Default = F.
                                
    --seed <int>          Set random seed.
                          Default = None.
            
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
    n_processes = 1
    output_dir = '.'
    file_root = 'super_percolator'
    static_mods = {'C': 57.02146}
    overwrite = False
    seed = None
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

    #df_all1['freq'] = df_all1['freq'].rank()

    train_power, std_power, true_power, discoveries = spf.do_iterative_svm_cv(df_all, folds=folds, Cs=[
        0.1, 1, 10], total_iter=total_iter, kernel=kernel, alpha=FDR_threshold, train_alpha=train_FDR_threshold, degree=degree, remove=remove)
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

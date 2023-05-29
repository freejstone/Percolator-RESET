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
import multiprocessing
import logging
import peptides
import sys
from tqdm import tqdm
import gzip
#########################################################################################################


def parse_static_mods(my_string):
  """
  Parse a static mod string (see USAGE) into a dictinoary.
  Key = amino acid, value = mass offset
  """
  return_value = {}

  for one_mod in my_string.split(","):
    words = one_mod.split(":")
    return_value[words[0]] = float(words[1])

  # Add in carbamidomethylation.
  if ("C" not in return_value):
    return_value["C"] = 57.02146

  return(return_value)
#########################################################################################################


def print_info(command_line, output_dir, file_root, overwrite, account_mods, search_file_narrow, search_file_open):
    #check if output directory exists, if not create and store log file there.
    if os.path.isdir(output_dir):
        if os.path.exists(path=output_dir + "/" + file_root + ".log.txt") and overwrite:
            os.remove(output_dir + "/" + file_root + ".log.txt")
            logging.basicConfig(filename=output_dir + "/" + file_root + ".log.txt",
                                level=logging.DEBUG, format='%(levelname)s: %(message)s')
        elif os.path.exists(path=output_dir + "/" + file_root + ".log.txt") and not overwrite:
            log_file = output_dir + "/" + file_root + ".log.txt"
            sys.exit("The file %s already exists and cannot be overwritten. Use --overwrite T to replace or choose a different output file name. \n" % (log_file))
        else:
            logging.basicConfig(filename=output_dir + "/" + file_root + ".log.txt",
                                level=logging.DEBUG, format='%(levelname)s: %(message)s')
    else:
        os.mkdir(output_dir)
        logging.basicConfig(filename=output_dir + "/" + file_root + ".log.txt",
                            level=logging.DEBUG, format='%(levelname)s: %(message)s')

    #print CPU info
    logging.info('CPU: ' + str(platform.platform()))
    sys.stderr.write('CPU: ' + str(platform.platform()) + " \n")

    #print version
    #logging.info('Version: ' + str(__version__))
    #sys.stderr.write('Version: ' + str(__version__) + " \n")

    #print date time info
    logging.info(str(datetime.datetime.now()))
    sys.stderr.write(str(datetime.datetime.now()) + " \n")

    #print command used
    logging.info('Command used: ' + command_line)
    sys.stderr.write('Command used: ' + command_line + "\n")

    sys.stderr.write("Successfully read in arguments. \n")
    logging.info("Successfully read in arguments")

    if not account_mods:
        logging.warning(
            "No longer accounting for variable modification. FDR control not guaranteed if variable modifications exist.")
        sys.stderr.write(
            "No longer accounting for variable modification. FDR control not guaranteed if variable modifications exist. \n")

    if type(search_file_narrow) == str:
        if (not os.path.isfile(search_file_narrow)):
            logging.info("The narrow search files does not exist.")
            sys.exit("The narrow search files does not exist. \n")
    if type(search_file_narrow) == list:
        if any(not os.path.isfile(i) for i in search_file_narrow):
            logging.info("One of the narrow search files does not exist.")
            sys.exit("One of the narrow search files does not exist. \n")
            
    if type(search_file_open) == str:
        if (not os.path.isfile(search_file_open)):
            logging.info("The open search files does not exist.")
            sys.exit("The open search files does not exist. \n")
    if type(search_file_open) == list:
        if any(not os.path.isfile(i) for i in search_file_open):
            logging.info("One of the open search files does not exist.")
            sys.exit("One of the open search files does not exist. \n")
    
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


def make_binidx_matchcount_map(mzs, fragment_min_mz, fragment_max_mz, bin_size):
    """
    Utility function for calc_proportion_fragments_incommon.
    Notionally bin the mzs from a list of mzs using the specified bin
    size (don't actually build the binned array), with specified lower
    and upper limits. Return a map from bin indexes to a count of
    fragments
    :param mzs: 
    :param fragment_min_mz:
    :param fragment_max_mz:
    :param bin_size:
    :return:  a map from bin index to a count of fragments in that bin
    """
    binidx_matchcount_map = {}
    for mz in mzs:
        if mz < fragment_min_mz or mz > fragment_max_mz:
            continue
        bin_idx = int((mz - fragment_min_mz) / bin_size)
        if bin_idx not in binidx_matchcount_map:
            binidx_matchcount_map[bin_idx] = 0
        binidx_matchcount_map[bin_idx] += 1
    return binidx_matchcount_map
###############################################################################


def calc_proportion_fragments_incommon(pepseq1, peptide1_mods, nterm1, cterm1,
                                       pepseq2, peptide2_mods, nterm2, cterm2,
                                       binsize, charge1, charge2):
    """
    Determine all the fragment mzs for each peptide. Bin the fragments from
    each peptide.  Calculate the fraction of fragments that fall into a bin
    with a fragment from the other peptide.

    :param pepseq1: First peptide sequence
    :param peptide1_mods: List of mass offsets, same length as peptide.
    :param nterm1: N-terminal modification mass of first peptide.
    :param cterm1: C-terminal modification mass of first peptide.
    :param pepseq2: Second peptide sequence
    :param peptide2_mods: List of mass offsets, same length as peptide.
    :param nterm2: N-terminal modification mass of second peptide.
    :param cterm2: C-terminal modification mass of second peptide.
    :param binsize: Size of m/z bins.
    :return: Fraction of matched peaks.
    """

    # Compute fragments and put them in bins.
    if charge1 in [1, 2]:
        fragment_charge1 = [1]
    elif charge1 >= 3:
        fragment_charge1 = [1, 2]

    if charge2 in [1, 2]:
        fragment_charge2 = [1]
    elif charge2 >= 3:
        fragment_charge2 = [1, 2]

    mzs_1 = peptides.calc_theoretical_peak_mzs(pepseq1, fragment_charge1, peptide1_mods,
                                               200, 3000, nterm1, cterm1)
    mzs_2 = peptides.calc_theoretical_peak_mzs(pepseq2, fragment_charge2, peptide2_mods,
                                               200, 3000, nterm2, cterm2)
    bin_count_map_1 = make_binidx_matchcount_map(mzs_1, 200, 3000, binsize)
    bin_count_map_2 = make_binidx_matchcount_map(mzs_2, 200, 3000, binsize)

    # Count matched bins.
    n_fragments_in_matched_bins = 0
    for binidx in bin_count_map_1:
        if binidx in bin_count_map_2:
            n_fragments_in_matched_bins += (bin_count_map_1[binidx]
                                            + bin_count_map_2[binidx])

    return float(n_fragments_in_matched_bins) / (len(mzs_1) + len(mzs_2))

###############################################################################


def parse_mods(pepseq_with_mods, static_mods):
    """
    Parse a modified peptide sequence string.

    :param pepseq_with_mods: Peptide string with interpolated bracketed mods.  
    :param static_mods: Dictionary of static mods. 
                        Key = amino acid, value = mass offset.
    :return: A list of amino acids, a list of modification values, and 
             the n-term and c-term terminal modification values.
    """
    aa_seq_list = []
    modifications = []
    nterm_delta = 0.0
    cterm_delta = 0.0

    aaseq_position = -1
    modpepseq_position = 0
    while modpepseq_position < len(pepseq_with_mods):
        my_char = pepseq_with_mods[modpepseq_position]
        if my_char.isalpha():

            # Create an unmodified amino acid and add it to the growing list.
            aa_seq_list.append(my_char)
            modifications.append(0.0)
            aaseq_position += 1
            modpepseq_position += 1
        elif (my_char == '['):
            end_pos = (pepseq_with_mods[modpepseq_position + 1:].index(']')
                       + modpepseq_position + 1)
            mod_mass = float(pepseq_with_mods[modpepseq_position + 1:end_pos])

            # Store a modification associated with the current position.
            # To handle when comet has a mass mod due to n- or c-term AND a variable mod on the same site. MSFragger/Tide do not appear to have this doubling up issue!
            if (modifications[aaseq_position] != 0.0):
                modifications[aaseq_position] += mod_mass
            else:
                modifications[aaseq_position] = mod_mass
            modpepseq_position = end_pos + 1
        else:
            sys.stderr.write("Invalid character (%s) at position %d.\n"
                             % (my_char, modpepseq_position))
            sys.exit(1)

    # Add in the static mods.
    for index in range(0, len(aa_seq_list)):
        amino = aa_seq_list[index]
        if (amino in static_mods):
          modifications[index] += static_mods[amino]
    if ("nterm" in static_mods):
      nterm_delta = static_mods["nterm"]
    if ("cterm" in static_mods):
      cterm_delta = static_mods["cterm"]

    return(aa_seq_list, modifications, nterm_delta, cterm_delta)
###############################################################################


def get_similarity(list1, charges1, list2, charges2, frag_bin_size=0.05, static_mods={'C': 57.02146}):

  peptide_out_1 = []
  peptide_out_2 = []
  similarity_out = []
  start_index = 0

  # Compare each pair of peptides.
  num_pairs = 0
  num_printed = 0
  for index1 in range(0, len(list1)):
    peptide1 = list1[index1]
    (unmodified_peptide1, peptide1_mods, nterm1, cterm1) \
        = parse_mods(peptide1, static_mods)

    for index2 in range(start_index, len(list2)):
      peptide2 = list2[index2]
      num_pairs += 1

      # Don't bother if they are the same peptide.
      if (peptide1.replace("I", "L") == peptide2.replace("I", "L")):
        peptide_out_1.append(peptide1)
        peptide_out_2.append(peptide2)
        similarity_out.append(1.0)
      else:

        (unmodified_peptide2, peptide2_mods, nterm2, cterm2) \
            = parse_mods(peptide2, static_mods)

        charge1 = charges1[index1]
        charge2 = charges2[index2]

        similarity = calc_proportion_fragments_incommon(
            unmodified_peptide1, peptide1_mods, nterm1, cterm1,
            unmodified_peptide2, peptide2_mods, nterm2, cterm2,
            frag_bin_size, charge1, charge2)

        num_printed += 1

        peptide_out_1.append(peptide1)
        peptide_out_2.append(peptide2)
        similarity_out.append(similarity)
  return(peptide_out_1, peptide_out_2, similarity_out)

#########################################################################################################


def filter_scan(search_df, thresh=0.05, frag_bin_size=0.05, static_mods={'C': 57.02146}):
    '''
    Parameters
    ----------
    search_df : Pandas Dataframe
        Dataframe containing just one scan, and multiple PSMs.
    thresh : float, optional
        The similarity score threshold used to filter neighbouring PSMs. The default is 0.05.
    frag_bin_size : float, optional
        Size of m/z bins used to determine the shared number of b- and y-ions. The default is 0.05.
    static_mods : dic, optional
        Dictionary containing all static modifications. The default is {'C':57.02146}.

    Returns
    -------
    drop_scan : list
        A list of booleans indicating which PSMs should be filtered.

    '''

    search_df = search_df.reset_index(drop=True)
    n_scans = search_df.shape[0]

    peptide_1 = [search_df['Peptide'].loc[0]]
    charge_1 = [search_df['charge'].loc[0]]
    drop_scan = [False]*n_scans

    #if any(search_df['sequence'][1:].str.contains(peptide_1[0])):
    for top in range(1, n_scans):
        peptide_2 = [search_df['Peptide'].loc[top]]
        charge_2 = [search_df['charge'].loc[top]]
        results = get_similarity(
            peptide_1, charge_1, peptide_2, charge_2, frag_bin_size, static_mods)
        if any(sim >= thresh for sim in results[2]):
            drop_scan[top] = True
        else:
            peptide_1 += peptide_2
            charge_1 += charge_2

    return(drop_scan)

###############################################################################


# wrapper is used to update the manager queue every time the wrapper is called
def wrapper(df, q, thresh, frag_bin_size, static_mods):
    '''
    a wrapper for filter_scan that can be used for multiprocessing
    '''
    result = filter_scan(df, thresh, frag_bin_size, static_mods)
    q.put(1)
    return(result)
###############################################################################


def filter_scan_subset(df, q, task_number, return_dict, tide_used, thresh, static_mods):
    '''
    Effectively calls a wrapper for filter_scan that can be used for multiprocessing
    '''
    if tide_used == 'tide':
        sys.stderr.write("Starting Process " + str(task_number) + " \n")
        logging.info("Starting Process " + str(task_number))
        results = df.groupby(['file', 'scan', "charge", "ExpMass"]).apply(
            lambda x: wrapper(x, q, thresh, 0.05, static_mods))
        results = results.sort_index(ascending=False)
        results = results.apply(pd.Series).stack().reset_index()
        return_dict[task_number] = results[0]

###############################################################################


def listener(q, list_of_df, tide_used):
    '''
    constantly checks to see if the queue q has been updated
    '''
    if tide_used == 'tide':
        pbar = tqdm(total=sum(
            [len(j.groupby(['file', 'scan', "charge", "ExpMass"])) for j in list_of_df]))
        for item in iter(q.get, None):
            pbar.update(item)

###############################################################################


def filter_narrow_open(narrow_target_decoys, open_target_decoys, thresh=0.05, n_processes=1, neighbour_remove=True, tide_used='tide', static_mods={'C': 57.02146}):
    '''
    Parameters
    ----------
    narrow_target_decoys : Pandas Dataframe
        Top 1 PSM for concatenated search of target-decoy database
    open_target_decoys : Pandas Dataframe
        Concatenated serch of target-decoy database
    open_top : int
        The number of top PSMs for each scan in the open search
    to be considered. The default is 2.
    thresh : float
        The threshold used to determine neighbouring peptides.
        The default is 0.05.
    n_processes : int
        Number of proccesses to be used. The default is 1.
    neighbour_remove : bool
        Whether we remove neighbours or not. The default is True.
    tide_used : bool
        Whether the input dataframes come from tide or not.
        The default is True.
    static_mods : dict
        Static modifications. The default is {'C':57.02146}.

    Returns
    -------
    target_decoys_all : Pandas Dataframe
        Combined PSMs from both narrow and open search, with neighbours
        for each scan number.
    '''
    open_target_decoys['n_o'] = 0
    narrow_target_decoys['n_o'] = 1

    open_target_decoys['rank'] = open_target_decoys['SpecId'].apply(
        lambda x: int(x[-1]))  # getting the rank
    narrow_target_decoys['rank'] = narrow_target_decoys['SpecId'].apply(
        lambda x: int(x[-1]))  # getting the rank
    narrow_target_decoys = narrow_target_decoys.loc[narrow_target_decoys['rank'] == 1, :]

    target_decoys_all = pd.concat([narrow_target_decoys, open_target_decoys]).reset_index(
        drop=True)  # combine narrow and open PSMs

    SpecId_split = target_decoys_all.SpecId.apply(
        lambda x: x.split('_'))  # split SpecId

    target_decoys_all[['target_decoy', 'file', 'scan', 'charge', 'rank']] = pd.DataFrame(
        SpecId_split.tolist())  # create new columns from SpecID

    target_decoys_all['charge'] = target_decoys_all['charge'].astype(int)

    target_decoys_all['Peptide'] = target_decoys_all['Peptide'].apply(
        lambda x: x[2:(len(x) - 2)])
    
    target_decoys_all['rank'] = target_decoys_all['rank'].astype(int)

    #drop duplicate PSMs that are subsequently found in the open search
    if tide_used == 'tide':
        target_decoys_all = target_decoys_all.drop_duplicates(
            ['file', 'scan', 'charge', 'ExpMass', 'Peptide'])

    #makes sure the PSMs from the narrow and open search are ordered by scan first, then by their score
    #indeed it makes sense to use xcorr_score over tailor score, as tailor_score uses candidate PSMs to
    #normalise the xcorr_score - this differs from narrow to open.
    if tide_used == 'tide':
        target_decoys_all['XCorr'] = round(target_decoys_all['XCorr'], 6)
        target_decoys_all = target_decoys_all.sort_values(
            by=['file', 'scan', 'charge', 'ExpMass', 'XCorr'], ascending=False)

    if not neighbour_remove:
        logging.info("Not removing neighbours.")
        sys.stderr.write("Not removing neighbours.\n")
        return(target_decoys_all)

    logging.info("Filtering for neighbours.")
    sys.stderr.write("Filtering for neighbours.\n")

    target_decoys_all = target_decoys_all[~(
        target_decoys_all.Peptide.str.contains('B|J|O|U|X|Z', regex=True))]
    target_decoys_all.reset_index(drop=True, inplace=True)

    if n_processes == 1:
        tqdm.pandas()
        if tide_used == 'tide':
            results = target_decoys_all.groupby(['file', 'scan', 'charge', 'ExpMass']).progress_apply(
                lambda x: filter_scan(x, thresh, 0.05, static_mods))  # apply filtering by experimental scan

        results = results.sort_index(ascending=False)
        results = results.apply(pd.Series).stack().reset_index()

        # create drop_scan column indicating which PSMs are being kept
        target_decoys_all['drop_scan'] = results[0]
    else:
        #if more than 1 thread, we partition the dataframe
        if tide_used == "tide":
            target_decoys_all['split_col'] = pd.qcut(
                target_decoys_all['scan'], n_processes)

        target_decoys_grouped = target_decoys_all.groupby(
            target_decoys_all.split_col)
        list_of_df = [0]*n_processes  # create a list of partitioned dataframes
        for i in range(len(target_decoys_all['split_col'].unique())):
            list_of_df[i] = target_decoys_grouped.get_group(
                target_decoys_all['split_col'].unique()[i])
            list_of_df[i].reset_index(drop=True, inplace=True)

        manager = multiprocessing.Manager()
        q = manager.Queue()  # creating instance of manager queue
        # creating manager dict variable that can be used to store changes to the argument by ALL processes at the same time
        return_dict = manager.dict()
        proc = multiprocessing.Process(
            target=listener, args=(q, list_of_df, tide_used))
        proc.start()  # start listening for updates to manager queue
        workers = [multiprocessing.Process(target=filter_scan_subset, args=(
            list_of_df[i], q, i, return_dict, tide_used, thresh, static_mods)) for i in range(n_processes)]  # now run each of the processes
        for worker in workers:
            worker.start()
        for worker in workers:
            worker.join()
        q.put(None)
        proc.join()

        for i in range(len(target_decoys_all['split_col'].unique())):
            target_decoys_all.loc[target_decoys_all['split_col'] == target_decoys_all['split_col'].unique(
            )[i], 'drop_scan'] = list(return_dict[i])  # update the main dataframe with the results from each process

    # drop the neighbours
    target_decoys_all = target_decoys_all[target_decoys_all['drop_scan'] == False]
    target_decoys_all.reset_index(drop=True, inplace=True)

    target_decoys_all.drop(
        ['target_decoy', 'file', 'scan', 'charge', 'drop_scan'], axis=1, inplace=True)

    return(target_decoys_all)

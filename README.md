# Percolator with FDR control

This contains the code for the associated manuscript, titled "Splitting the decoys into training and estimating: how to train a Percolator-like post-processor and maintain control of the false discovery rate". 

The important folder is *functions* which contains the actual scripts that implements our procedure using Python, from the command line. The main python script is the file, `do_FDR_percolator.py`, and the auxiliary functions can be found under `percolator_functions.py` and `utility_functions.py`. This python implementation may be used as interim until the development in C (Percolator). For further details on implementation, see the USAGE under `do_FDR_percolator.py`.

A note to Lukas: for the development in C, please refer to the link [here](https://unisydneyedu-my.sharepoint.com/:f:/g/personal/jfre0619_uni_sydney_edu_au/Et239lX9DY9LjFA4GtfyHvQBMDRZbCc8dcU2QDaE_2mDPw?e=gUJDXR) which contains test files and their associated results using the current Python implementation.

## Description of files

* run_tide_index.sh: Digests protein database and creates 10 random decoy databases for the PRIDE-20 dataset.
* narrow_search.sh: Conducts 10 concatenated narrow searches for each random decoy database, for each PRIDE-20 dataset.
* open_search.sh: Conducts 10 concatenated open searches for each random decoy database, for each PRIDE-20 dataset.
* open_search_sep.sh: Conducts 10 separate open searches for each random decoy database, for each PRIDE-20 dataset.
* narrow_search_sep.sh: Conducts 10 separate narrow searches for each random decoy database, for each PRIDE-20 dataset.
* make_pin.sh: Produces pin files for open_search.sh and narrow_search.sh.
* enzint_top1.R: Selects the top 1 PSM and removes the EnzInt feature for subsequent use by Percolator.
* run_narrow_percolator.sh: Conducts Percolator on search files from narrow_search.sh.
* run_open_percolator.sh: Conducts Percolator on search files from open_search.sh.
* run_single_decoy.sh: Conducts Percolator-RESET files from on files from narrow_search.sh and open_search.sh.
* run_extra_decoy.sh: Conducts Percolator-RESET+ files from on files from narrow_search_sep.sh and open_search_sep.sh.
* run_single_decoy_var.sh: Conducts Percolator-RESET files from on files from narrow_search.sh and open_search.sh, with varying seeds.
* run_extra_decoy_var.sh: Conducts Percolator-RESET+ files from on files from narrow_search_sep.sh and open_search_sep.sh, with varying seeds.
* percolator_summary.R: Does peptide/stem-level competition for Percolator.

For files under entrapment_runs/:
* combine_fasta.R: Prepares the in-sample and entrapment sequences with decreasing in-sample proportion (100\%, 75\%, 50\%, 25\%).
* tide-index-trypsin.sh: Runs Tide-index to create 100 target-decoy databases using each of the generated fasta files from combine_fasta.R
* tide-search.sh: Runs concatenated narrow search for each of the databases given by tide-index-trypsin.sh.
* tide-search-sep.sh: Runs separate narrow search for each of the databases given by tide-index-trypsin.sh.
* tide-search-open.sh: Runs concatenated open search for each of the databases given by tide-index-trypsin.sh.
* tide-search-sep-open.sh: Runs separate open search for each of the databases given by tide-index-trypsin.sh.
* fix_enzit_tide.R: Selects the top 1 PSM and removes the EnzInt feature for subsequent use by Percolator.
* percolator.sh: Conducts Percolator on search files from tide-search.sh and tide-search-open.sh.
* run_single_decoy.sh: Conducts Percolator-RESET files from on files from tide-search.sh and tide-search-open.sh.
* run_extra_decoy.sh: Conducts Percolator-RESET+ files from on files from tide-search.sh and tide-search-open.sh.
* percolator_results.R: Does peptide/stem-level competition for Percolator.
* percolator_with_fdr_results.R: Summarises results from our method Percolator RESET.

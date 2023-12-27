## Description of files under code_for_paper/

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
* create_decoy_database.R: Produces the target-decoy protein fasta files (1 reversed and 1 cyclically permuted) used with MSFragger searches.
* get_frag_tol.R: Adds the fragment tolerances determined by Tide to the parameter files used with MSFragger.
* narrow_search_msfragger.sh: Do narrow searches using MSFragger.
* open_search_msfragger.sh: Do open searches using MSFragger.
* remove_duplicate_decoys.R: Removing peptides that are both 'targets' and 'decoys' and are thus ill-defined from MSFragger searches.
* run_percolator_msfragger.sh: Running Percolator on the MSFragger searches.
* run_single_msfragger.sh: Running Percolator-RESET on MSFragger searches.
* run_extra_decoy_msfragger.sh: Running Percolator-RESET+ on MSFragger searches.
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
* run_single_decoy_no_pair.sh: Conducts Percolator-RESET files from on files from tide-search.sh and tide-search-open.sh without target-decoy pairing.
* run_extra_decoy_no_pair.sh: Conducts Percolator-RESET+ files from on files from tide-search.sh and tide-search-open.sh without target-decoy pairing.
* percolator_with_fdr_results.R: Summarises results from our method Percolator RESET.
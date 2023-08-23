#!/bin/bash

#varies the training decoy set

# Specify the file path
file_path="./datasets"

folders=($(ls -d "$file_path"/*/))
values=($(seq 0 19))

for value in "${values[@]}"; do
    for folder in "${folders[@]}"; do
        python functions/do_FDR_percolator.py --seed "${value}" --FDR_threshold 0.3 --overwrite T --output_dir "${folder}crux-output" --file_root "open_extra_decoy_seed_${value}_ind" "${folder}crux-output/open_50_sep.tide-search.target.pin,${folder}crux-output/open_51_sep.tide-search.target.pin,${folder}crux-output/open_50_sep.tide-search.decoy.pin,${folder}crux-output/open_51_sep.tide-search.decoy.pin" "${folder}index-0/tide-index.peptides.txt,${folder}index-1/tide-index.peptides.txt"
    done
done

for value in "${values[@]}"; do
    for folder in "${folders[@]}"; do
        python functions/do_FDR_percolator.py --seed "${value}" --FDR_threshold 0.3 --overwrite T --output_dir "${folder}crux-output" --file_root "narrow_extra_decoy_seed_${value}_ind" "${folder}crux-output/narrow_50_sep.tide-search.target.pin,${folder}crux-output/narrow_51_sep.tide-search.target.pin,${folder}crux-output/narrow_50_sep.tide-search.decoy.pin,${folder}crux-output/narrow_51_sep.tide-search.decoy.pin" "${folder}index-0/tide-index.peptides.txt,${folder}index-1/tide-index.peptides.txt"
    done
done


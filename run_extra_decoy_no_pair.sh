#!/bin/bash

# Specify the file path
file_path="../datasets"

folders=($(ls -d "$file_path"/*/))
values=($(seq 0 9))

for value in "${values[@]}"; do
    for folder in "${folders[@]}"; do
        value1=$(( (value + 1) % 10 ))
        python do_FDR_percolator.py --seed "${value}" --FDR_threshold 0.3 --overwrite T --narrow F --pair F --output_dir "${folder}crux-output" --file_root "open_extra_decoy_${value}_no_pair" "${folder}crux-output/open_dups_removed_${value}_sep.tide-search.target.pin,${folder}crux-output/open_dups_removed_${value1}_sep.tide-search.target.pin,${folder}crux-output/open_dups_removed_${value}_sep.tide-search.decoy.pin,${folder}crux-output/open_dups_removed_${value1}_sep.tide-search.decoy.pin"
    done
done

for value in "${values[@]}"; do
    for folder in "${folders[@]}"; do
        value1=$(( (value + 1) % 10 ))
        python do_FDR_percolator.py --seed "${value}" --FDR_threshold 0.3 --overwrite T --narrow T --pair F --output_dir "${folder}crux-output" --file_root "narrow_extra_decoy_${value}_no_pair" "${folder}crux-output/narrow_dups_removed_${value}_sep.tide-search.target.pin,${folder}crux-output/narrow_dups_removed_${value1}_sep.tide-search.target.pin,${folder}crux-output/narrow_dups_removed_${value}_sep.tide-search.decoy.pin,${folder}crux-output/narrow_dups_removed_${value1}_sep.tide-search.decoy.pin"
    done
done


#!/bin/bash

# Specify the file path
file_path="./datasets"

folders=($(ls -d "$file_path"/*/))
values=($(seq 0 9))

for value in "${values[@]}"; do
    for folder in "${folders[@]}"; do
        value1=$(( (value + 1) % 10 ))
        python3 -m percolator_RESET --seed "$((${value}+19122023))" --FDR_threshold 0.3 --overwrite T --initial_dir TailorScore --score TailorScore --pair F --output_dir "${folder}crux-output" --file_root "open_extra_decoy_${value}_no_pair_fix" "${folder}crux-output/open_dups_removed_${value}_sep.tide-search.target.pin,${folder}crux-output/open_dups_removed_${value1}_sep.tide-search.target.pin,${folder}crux-output/open_dups_removed_${value}_sep.tide-search.decoy.pin,${folder}crux-output/open_dups_removed_${value1}_sep.tide-search.decoy.pin"
    done
done

for value in "${values[@]}"; do
    for folder in "${folders[@]}"; do
        value1=$(( (value + 1) % 10 ))
        python3 -m percolator_RESET --seed "$((${value}+19122023))" --FDR_threshold 0.3 --overwrite T --initial_dir TailorScore --score TailorScore --pair F --output_dir "${folder}crux-output" --file_root "narrow_extra_decoy_${value}_no_pair_fix" "${folder}crux-output/narrow_dups_removed_${value}_sep.tide-search.target.pin,${folder}crux-output/narrow_dups_removed_${value1}_sep.tide-search.target.pin,${folder}crux-output/narrow_dups_removed_${value}_sep.tide-search.decoy.pin,${folder}crux-output/narrow_dups_removed_${value1}_sep.tide-search.decoy.pin"
    done
done


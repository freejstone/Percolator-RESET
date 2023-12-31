#!/bin/bash

# Specify the file path
file_path="./datasets"

folders=($(ls -d "$file_path"/*/))
values=($(seq 0 1))
runs=($(seq 0 4))
its=($(seq 0 19))
initial=15112023


for it in "${its[@]}"; do
    count=0
    for value in "${values[@]}"; do
        folder=${folders[it]}
        for run in "${runs[@]}"; do
            initial=$((initial + 1))
            python3 -m percolator_RESET --seed "${initial}" --FDR_threshold 0.05 --overwrite T --pair F --initial_dir hyperscore --output_dir "${folder}msfragger_output" --file_root "narrow_single_decoy_${count}_ind_fix" "${folder}msfragger_output/narrow_dups_removed_${value}.pin"
            count=$((count + 1))
        done
    done
done

for it in "${its[@]}"; do
    count=0
    for value in "${values[@]}"; do
        folder=${folders[it]} 
        for run in "${runs[@]}"; do
            initial=$((initial + 1))
            python3 -m percolator_RESET --seed "${initial}" --FDR_threshold 0.05 --overwrite T --pair F --initial_dir hyperscore --output_dir "${folder}msfragger_output" --file_root "open_single_decoy_${count}_ind_fix" "${folder}msfragger_output/open_dups_removed_${value}.pin"
            count=$((count + 1))
        done
    done
done





#!/bin/bash

# Specify the file path
file_path="./datasets"

folders=($(ls -d "$file_path"/*/))
values=($(seq 0 0))
runs=($(seq 0 9))
its=($(seq 0 19))
initial=15112023


for it in "${its[@]}"; do
    count=0
    for value in "${values[@]}"; do
        folder=${folders[it]}
        for run in "${runs[@]}"; do
            initial=$((initial + 1))
            python3 -m percolator_RESET --seed "${initial}" --FDR_threshold 0.05 --mult 2 --overwrite T --pair F --initial_dir hyperscore --output_dir "${folder}msfragger_output" --file_root "narrow_extra_decoy_${count}_ind_fix" "${folder}msfragger_output/narrow_dups_removed.pin"
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
            python3 -m percolator_RESET --seed "${initial}" --FDR_threshold 0.05 --mult 2 --overwrite T --pair F --initial_dir hyperscore --output_dir "${folder}msfragger_output" --file_root "open_extra_decoy_${count}_ind_fix" "${folder}msfragger_output/open_dups_removed.pin"
            count=$((count + 1))
        done
    done
done





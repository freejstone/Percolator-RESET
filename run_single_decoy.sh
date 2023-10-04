#!/bin/bash

# Specify the file path
file_path="../datasets"

folders=($(ls -d "$file_path"/*/))
values=($(seq 0 9))
its=($(seq 0 19))

for it in "${its[@]}"; do
    for value in "${values[@]}"; do
        folder=${folders[it]} 
        python do_FDR_percolator.py --seed "${value}" --FDR_threshold 0.3 --overwrite T --narrow T --output_dir "${folder}crux-output" --file_root "narrow_single_decoy_${value}_ind" "${folder}crux-output/narrow_5${value}.make-pin.pin" "${folder}index-${value}/tide-index.peptides.txt"
    done
done


for it in "${its[@]}"; do
    for value in "${values[@]}"; do
        folder=${folders[it]} 
        python do_FDR_percolator.py --seed "${value}" --FDR_threshold 0.3 --overwrite T --narrow F --output_dir "${folder}crux-output" --file_root "open_single_decoy_${value}_ind" "${folder}crux-output/open_5${value}.make-pin.pin" "${folder}index-${value}/tide-index.peptides.txt"
    done
done





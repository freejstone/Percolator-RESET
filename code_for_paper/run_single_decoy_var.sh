#!/bin/bash

#varies the training decoy set

# Specify the file path
file_path="./datasets"

folders=($(ls -d "$file_path"/*/))
decoys=($(seq 0 1))
values=($(seq 0 9))
its=($(seq 0 19))

for it in "${its[@]}"; do
    for decoy in "${decoys[@]}"; do
        for value in "${values[@]}"; do
            folder=${folders[it]} 
            python -m percolator_RESET --seed "${value}" --FDR_threshold 0.3 --overwrite T --initial_dir TailorScore --score TailorScore --output_dir "${folder}crux-output" --file_root "narrow_single_decoy_${decoy}_seed_${value}_ind_fix" "${folder}crux-output/narrow_5${decoy}.make-pin.pin" "${folder}index-${decoy}/tide-index.peptides.txt"
        done
    done
done

for it in "${its[@]}"; do
    for decoy in "${decoys[@]}"; do
        for value in "${values[@]}"; do
            folder=${folders[it]} 
            python -m percolator_RESET --seed "${value}" --FDR_threshold 0.3 --overwrite T --initial_dir TailorScore --score TailorScore --output_dir "${folder}crux-output" --file_root "open_single_decoy_${decoy}_seed_${value}_ind_fix" "${folder}crux-output/open_5${decoy}.make-pin.pin" "${folder}index-${decoy}/tide-index.peptides.txt"
        done
    done
done




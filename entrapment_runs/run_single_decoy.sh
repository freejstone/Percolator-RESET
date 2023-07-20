#!/bin/bash

total=100

section1() {
    sizes=('025')
    values=($(seq 0 99))
    repeats=($(seq 0 3))


    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do 
                python ../functions/do_FDR_percolator.py --FDR_threshold 0.3 --overwrite T --output_dir percolator_with_fdr_results --file_root "narrow_single_decoy_${size}_${counter}" "tide-search/${size}/dcy_${size}_${value}.tide-search.pin" "tide-index/index-${size}-${value}/tide-index.peptides.txt"
                counter=$((counter + 1))
            done
        done
    done

    for size in "${sizes[@]}"; do
        counter=0
       for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do 
                python ../functions/do_FDR_percolator.py --FDR_threshold 0.3 --overwrite T --output_dir percolator_with_fdr_results --file_root "open_single_decoy_${size}_${counter}" "tide-search/${size}/open_top5_${size}_${value}.tide-search.pin" "tide-index/index-${size}-${value}/tide-index.peptides.txt"
                counter=$((counter + 1))
            done
        done
    done
}


section2() {
    sizes=('half')
    values=($(seq 0 99))
    repeats=($(seq 0 3))


    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do 
                python ../functions/do_FDR_percolator.py --FDR_threshold 0.3 --overwrite T --output_dir percolator_with_fdr_results --file_root "narrow_single_decoy_${size}_${counter}" "tide-search/${size}/dcy_${size}_${value}.tide-search.pin" "tide-index/index-${size}-${value}/tide-index.peptides.txt"
                counter=$((counter + 1))
            done
        done
    done

    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do 
                python ../functions/do_FDR_percolator.py --FDR_threshold 0.3 --overwrite T --output_dir percolator_with_fdr_results --file_root "open_single_decoy_${size}_${counter}" "tide-search/${size}/open_top5_${size}_${value}.tide-search.pin" "tide-index/index-${size}-${value}/tide-index.peptides.txt"
                counter=$((counter + 1))
            done
        done
    done
}


section3() {
    sizes=('075')
    values=($(seq 0 99))
    repeats=($(seq 0 4))


    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do 
                python ../functions/do_FDR_percolator.py --FDR_threshold 0.3 --overwrite T --output_dir percolator_with_fdr_results --file_root "narrow_single_decoy_${size}_${counter}" "tide-search/${size}/dcy_${size}_${value}.tide-search.pin" "tide-index/index-${size}-${value}/tide-index.peptides.txt"
                counter=$((counter + 1))
            done
        done
    done

    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do 
                python ../functions/do_FDR_percolator.py --FDR_threshold 0.3 --overwrite T --output_dir percolator_with_fdr_results --file_root "open_single_decoy_${size}_${counter}" "tide-search/${size}/open_top5_${size}_${value}.tide-search.pin" "tide-index/index-${size}-${value}/tide-index.peptides.txt"
                counter=$((counter + 1))
            done
        done
    done
}

section4() {
    sizes=('full')
    values=($(seq 0 99))
    repeats=($(seq 0 4))


    for size in "${sizes[@]}"; do
        counter=400
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do 
                python ../functions/do_FDR_percolator.py --FDR_threshold 0.3 --overwrite T --output_dir percolator_with_fdr_results --file_root "narrow_single_decoy_${size}_${counter}" "tide-search/${size}/dcy_${size}_${value}.tide-search.pin" "tide-index/index-${size}-${value}/tide-index.peptides.txt"
                counter=$((counter + 1))
            done
        done
    done

    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do 
                python ../functions/do_FDR_percolator.py --FDR_threshold 0.3 --overwrite T --output_dir percolator_with_fdr_results --file_root "open_single_decoy_${size}_${counter}" "tide-search/${size}/open_top5_${size}_${value}.tide-search.pin" "tide-index/index-${size}-${value}/tide-index.peptides.txt"
                counter=$((counter + 1))
            done
        done
    done
}

# Check command line arguments and execute the corresponding section

if [[ $# -eq 0 ]]; then
    echo "No section specified. Usage: ./run_single_decoy.sh [section]"
    exit 1
fi

section=$1

case $section in
    section1)
        section1
        ;;
    section2)
        section2
        ;;
    section3)
        section3
        ;;
    section4)
        section4
        ;;
    *)
        echo "Invalid section. Available sections: section1, section2, section3, section4"
        exit 1
        ;;
esac

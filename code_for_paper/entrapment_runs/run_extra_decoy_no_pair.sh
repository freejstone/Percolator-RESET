#!/bin/bash

total=100

section1() {
    sizes=('025')
    values=($(seq 0 99))
    repeats=($(seq 0 1))


    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do
                value1=$(( (value + 1) % total ))
                python3 -m percolator_RESET --FDR_threshold 0.3 --overwrite T --pair F --mult 2 --initial_dir TailorScore --score TailorScore --seed "${counter}" --output_dir percolator_with_fdr_results --file_root "narrow_extra_decoy_${size}_${counter}_no_pair_fix" "tide-search/${size}/dcy_sep_${size}_${value}.tide-search.target.pin,tide-search/${size}/dcy_sep_${size}_${value1}.tide-search.target.pin,tide-search/${size}/dcy_sep_${size}_${value}.tide-search.decoy.pin,tide-search/${size}/dcy_sep_${size}_${value1}.tide-search.decoy.pin"
                counter=$((counter + 1))
            done
        done
    done

    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do
                value1=$(( (value + 1) % total ))
                python3 -m percolator_RESET --FDR_threshold 0.3 --overwrite T --pair F --mult 2 --initial_dir TailorScore --score TailorScore --seed "${counter}" --output_dir percolator_with_fdr_results --file_root "open_extra_decoy_${size}_${counter}_no_pair_fix" "tide-search/${size}/open_top5_sep_${size}_${value}.tide-search.target.pin,tide-search/${size}/open_top5_sep_${size}_${value1}.tide-search.target.pin,tide-search/${size}/open_top5_sep_${size}_${value}.tide-search.decoy.pin,tide-search/${size}/open_top5_sep_${size}_${value1}.tide-search.decoy.pin"
                counter=$((counter + 1))
            done
        done
    done
}


section2() {
    sizes=('half')
    values=($(seq 0 99))
    repeats=($(seq 0 1))


    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do
                value1=$(( (value + 1) % total ))
                python3 -m percolator_RESET --FDR_threshold 0.3 --overwrite T --pair F --mult 2 --initial_dir TailorScore --score TailorScore --seed "${counter}" --output_dir percolator_with_fdr_results --file_root "narrow_extra_decoy_${size}_${counter}_no_pair_fix" "tide-search/${size}/dcy_sep_${size}_${value}.tide-search.target.pin,tide-search/${size}/dcy_sep_${size}_${value1}.tide-search.target.pin,tide-search/${size}/dcy_sep_${size}_${value}.tide-search.decoy.pin,tide-search/${size}/dcy_sep_${size}_${value1}.tide-search.decoy.pin" 
                counter=$((counter + 1))
            done
        done
    done

    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do
                value1=$(( (value + 1) % total ))
                python3 -m percolator_RESET --FDR_threshold 0.3 --overwrite T --pair F --mult 2 --initial_dir TailorScore --score TailorScore --seed "${counter}" --output_dir percolator_with_fdr_results --file_root "open_extra_decoy_${size}_${counter}_no_pair_fix" "tide-search/${size}/open_top5_sep_${size}_${value}.tide-search.target.pin,tide-search/${size}/open_top5_sep_${size}_${value1}.tide-search.target.pin,tide-search/${size}/open_top5_sep_${size}_${value}.tide-search.decoy.pin,tide-search/${size}/open_top5_sep_${size}_${value1}.tide-search.decoy.pin"
                counter=$((counter + 1))
            done
        done
    done
}


section3() {
    sizes=('075')
    values=($(seq 0 99))
    repeats=($(seq 0 1))


    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do
                value1=$(( (value + 1) % total ))
                python3 -m percolator_RESET --FDR_threshold 0.3 --overwrite T --pair F --mult 2 --initial_dir TailorScore --score TailorScore --seed "${counter}" --output_dir percolator_with_fdr_results --file_root "narrow_extra_decoy_${size}_${counter}_no_pair_fix" "tide-search/${size}/dcy_sep_${size}_${value}.tide-search.target.pin,tide-search/${size}/dcy_sep_${size}_${value1}.tide-search.target.pin,tide-search/${size}/dcy_sep_${size}_${value}.tide-search.decoy.pin,tide-search/${size}/dcy_sep_${size}_${value1}.tide-search.decoy.pin"
                counter=$((counter + 1))
            done
        done
    done

    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do
                value1=$(( (value + 1) % total ))
                python3 -m percolator_RESET --FDR_threshold 0.3 --overwrite T --pair F --mult 2 --initial_dir TailorScore --score TailorScore --seed "${counter}" --output_dir percolator_with_fdr_results --file_root "open_extra_decoy_${size}_${counter}_no_pair_fix" "tide-search/${size}/open_top5_sep_${size}_${value}.tide-search.target.pin,tide-search/${size}/open_top5_sep_${size}_${value1}.tide-search.target.pin,tide-search/${size}/open_top5_sep_${size}_${value}.tide-search.decoy.pin,tide-search/${size}/open_top5_sep_${size}_${value1}.tide-search.decoy.pin"
                counter=$((counter + 1))
            done
        done
    done
}

section4() {
    sizes=('full')
    values=($(seq 0 99))
    repeats=($(seq 0 1))

    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do
                value1=$(( (value + 1) % total ))
                python3 -m percolator_RESET --FDR_threshold 0.3 --overwrite T --pair F --mult 2 --initial_dir TailorScore --score TailorScore --seed "${counter}" --output_dir percolator_with_fdr_results --file_root "open_extra_decoy_${size}_${counter}_no_pair_fix" "tide-search/${size}/open_top5_sep_${size}_${value}.tide-search.target.pin,tide-search/${size}/open_top5_sep_${size}_${value1}.tide-search.target.pin,tide-search/${size}/open_top5_sep_${size}_${value}.tide-search.decoy.pin,tide-search/${size}/open_top5_sep_${size}_${value1}.tide-search.decoy.pin"
                counter=$((counter + 1))
            done
        done
    done

    for size in "${sizes[@]}"; do
        counter=0
        for repeat in "${repeats[@]}"; do
            for value in "${values[@]}"; do
                value1=$(( (value + 1) % total ))
                python3 -m percolator_RESET --FDR_threshold 0.3 --overwrite T --pair F --mult 2 --initial_dir TailorScore --score TailorScore --seed "${counter}" --output_dir percolator_with_fdr_results --file_root "narrow_extra_decoy_${size}_${counter}_no_pair_fix" "tide-search/${size}/dcy_sep_${size}_${value}.tide-search.target.pin,tide-search/${size}/dcy_sep_${size}_${value1}.tide-search.target.pin,tide-search/${size}/dcy_sep_${size}_${value}.tide-search.decoy.pin,tide-search/${size}/dcy_sep_${size}_${value1}.tide-search.decoy.pin" 
                counter=$((counter + 1))
            done
        done
    done
}


# Check command line arguments and execute the corresponding section

if [[ $# -eq 0 ]]; then
    echo "No section specified. Usage: ./run_extra_decoy_no_pair.sh [section]"
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
    section5)
        section5
        ;;
    *)
        echo "Invalid section. Available sections: section1, section2, section3, section4"
        exit 1
        ;;
esac
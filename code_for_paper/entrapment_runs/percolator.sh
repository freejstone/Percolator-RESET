#!/bin/bash

sizes=('025' 'half' '075' 'full')
values=($(seq 0 99))

for size in "${sizes[@]}"; do
    for value in "${values[@]}"; do 
        crux percolator --only-psms T --overwrite T --fileroot open_top1_enzint_zero_${size}_${value} --tdc F --output-dir percolator_results tide-search/${size}/open_top1_enzint_zero_${size}_${value}.tide-search.pin
    done
done

for size in "${sizes[@]}"; do
    for value in "${values[@]}"; do 
        crux percolator --only-psms T --overwrite T --fileroot dcy_enzint_zero_${size}_${value} --tdc F --output-dir percolator_results tide-search/${size}/dcy_enzint_zero_${size}_${value}.tide-search.pin
    done
done
#!/bin/bash

declare base_dir="$(dirname $(dirname $(realpath $0)))"
cd ${base_dir}

mkdir -p build
mkdir -p results

cd ${base_dir}/build
set -e
cmake ..
make -j
set +e

n_list=(10000000)
./benchmark_peak_space "${n_list[@]}"

cd ${base_dir}/scripts
for n in "${n_list[@]}"; do
    python3 plot_space_results.py ../results/peak_space_$n.csv ../results/peak_space.pdf
done

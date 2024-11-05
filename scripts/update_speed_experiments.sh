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

n_list=(100000000)
./benchmark_update_speed "${n_list[@]}"

cd ${base_dir}/scripts
for n in "${n_list[@]}"; do
    python3 plot_results.py ../results/update_speed_$n.csv ../results/update_speed_$n.pdf
done

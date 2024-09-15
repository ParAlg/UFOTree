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

n_list=(1000 10000 100000)
./benchmark_query_speed "${n_list[@]}"

cd ${base_dir}/scripts
for n in "${n_list[@]}"; do
    python3 plot_results.py ../results/query_speed_$n.csv ../results/query_speed_$n.pdf
done

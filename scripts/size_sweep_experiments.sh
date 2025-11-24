#!/bin/bash

declare base_dir="$(dirname $(dirname $(realpath $0)))"
cd ${base_dir}

mkdir -p build
mkdir -p results

cd ${base_dir}/build
set -e
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j
set +e

./parallel_benchmark_n_sweep

cd ${base_dir}/scripts
python3 plot_n_sweep_parallel_updates.py

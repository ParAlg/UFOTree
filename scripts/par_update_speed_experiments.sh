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

export PARLAY_NUM_THREADS=192

N=10000000
K=1000000
./parallel_benchmark "$N" "$K"
./parallel_graph_benchmark "$K"

cd ${base_dir}/scripts
python3 plot_par_update_results.py ../results/parallel_update_speed_${N}_${K}.csv ../results/parallel_update_speed_graph.csv ../results/par_update_speed.pdf

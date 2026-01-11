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

N=10000000
./benchmark_update_speed "$N"
./benchmark_graph_update_speed

cd ${base_dir}/scripts
python3 plot_seq_update_results.py ../results/update_speed_$N.csv ../results/update_speed_graph.csv ../results/update_speed.pdf
# python3 plot_results.py ../results/update_speed_$N.csv ../results/update_speed_$N.csv ../results/update_speed.pdf

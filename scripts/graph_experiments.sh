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

./benchmark_graph_update_speed
python3 plot_results.py ../results/update_speed_graph.csv ../results/update_speed_graph.pdf

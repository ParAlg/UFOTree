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
NQ=1000000
Q=1000000

./benchmark_updates_diameter_sweep "$N"
./benchmark_queries_diameter_sweep "$NQ" "$Q"

cd ${base_dir}/scripts
python3 plot_diameter_sweep_sequential.py ../results/diameter_sweep_update_$N.csv ../results/diameter_sweep_path_query_$NQ.csv ../results/diameter_sweep_conn_query_$NQ.csv ../results/diameter_sweep.pdf

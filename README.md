# UFO Trees: Practical and Provably-Efficient Parallel Batch-Dynamic Trees (PPoPP 2026)

## Getting Started Guide

### Prerequisites

You will need the following tools installed on your system before beginning the build process:

  * **CMake**: Version **3.15** or higher is required for configuration.
  * **G++**: This project is validated for **g++ version 11.4.0**.
  * **Python**: Version **3.9** or higher.
      * Required Python packages: `matplotlib`, `numpy`, `pandas`.

#### Verifying Your Environment

You can check your current versions by running the following commands in your terminal:

```bash
cmake --version
g++ --version
python3 --version
```

### Building the Artifact

Building the artifact can be done with the following commands in your terminal. Alternatively, you can just run one of the experimental scripts which will automatically create the necessary directories, build the project, run it, and plot the results.

```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```

### Minimum Example

To get started run the shell script:

```bash
scripts/mini_experiments.sh
```

This should take no more than a minute. It will attempt to build the project and run a small experiment on very small synthetic trees with $n=1000$. The results are written into `results/update_speed_1000.csv` and are plotted in `results/update_speed_mini.pdf`.
* _Note: This plot duplicates the results in the top and bottom. In our full experiments we will use the top row for synthetic data and the bottom row for real-world data._



## Step-by-Step Instructions

Here we give detailed instructions for how to reproduce the results for each of the main experiments of our paper. First we describe how to download and set up the graph datasets that we use. Then we describe each experiment in detail.
This table summarizes the scripts to reproduce each experiment:


| Experiment | Figure | Script (in `scripts/`) | Result PDF (in `results/`) |
| :--- | :--- | :--- | :--- |
| Sequential Update Speed | Figure 5 | `update_speed_experiments.sh` | `update_speed.pdf` |
| Diameter Sweep | Figure 6 | `diameter_sweep_experiments.sh` | `diameter_sweep.pdf` |
| Parallel Update Speed | Figure 7 | `par_update_speed_experiments.sh` | `par_update_speed.pdf` |
| Input Size Scaling | Figure 8 | `size_sweep_experiments.sh` | `n_sweep.pdf` |

### Graph Datasets

Download the graph datasets provided, and store them (or symlink them) in the directory `build/graphs/`. The graphs used in our experiments by default are listed as follows, and must have the same names for the artifact to work as is:

* RoadUSA_sym.bin
* enwiki_sym.bin
* stackoverflow_sym.bin
* twitter_sym.bin


### Sequential Update Speed Experiments (Figure 5)

This experiment is completely sequential and only requires a single thread (although some parallelism is used to generate the test data). Run the shell script:

```bash
scripts/update_speed_experiments.sh
```

* Results written to: `results/update_speed_10000000.csv` and `results/update_speed_graph.csv`.
* Results plotted in: `results/update_speed.pdf`.

Customize This Experiment:
* These experiments take the average time over $3$ trials for each case. You can change this value by modifying the $4$-th value in each tuple entry in the `test_cases` lists in `benchmark/dynamic_trees/benchmark_update_speed.cpp` and `benchmark/graph_algorithms/benchmark_graph_update_speed.cpp`. To save time we recommend setting this to $1$ for all cases.
* To change the value of $n$ for the synthetic trees in this experiment, you can modify the `n_list` variable in `scripts/update_speed_experiments.sh`.
* To change the real-world graphs that are used in this experiment, you can modify the `test_cases` variable in `benchmark/graph_algorithms/benchmark_graph_update_speed.cpp`.
* To change the data structures that are tested in these experiments you can comment out the relevant lines in `benchmark/dynamic_trees/benchmark_update_speed.cpp` and `benchmark/graph_algorithms/benchmark_graph_update_speed.cpp`. To save the most time we recommend not running the RC-Tree and Topology Tree baselines as they are significantly slower on some inputs.


### Diameter Sweep and Query Experiments (Figure 6)

This experiment is completely sequential and only requires a single thread (although some parallelism is used to generate the test data). Run the shell script:

```bash
scripts/diameter_sweep_experiments.sh
```

* Results written to: `results/diameter_sweep_update_10000000.csv`, `results/diameter_sweep_conn_query_1000000.csv`, and `results/diameter_sweep_path_query_1000000.csv`.
* Results plotted in: `results/diameter_sweep.pdf`.

Customize This Experiment:
* These experiments take the average time over $3$ trials for each case. You can change this value by modifying the $4$-th value in each tuple entry in the `test_cases` lists in `benchmark/dynamic_trees/benchmark_updates_diameter_sweep.cpp` and `benchmark/dynamic_trees/benchmark_queries_diameter_sweep.cpp`. To save time we recommend setting this to $1$ for all cases.
* You can change the values of `N`, `NQ`, and `Q` in `scripts/diameter_sweep_experiments.sh`.
  * `N` is the size of the trees for the update speed diameter sweep experiment.
  * `NQ` is the size of the trees for the query speed diameter sweep experiment.
  * `Q` is the number of queries in the query speed diameter sweep experiment.
* To change the values of Alpha in these diameter sweeps, you can modify the first entries in each tuple in the `test_cases` lists in both `benchmark/dynamic_trees/benchmark_updates_diameter_sweep.cpp` and `benchmark/dynamic_trees/benchmark_queries_diameter_sweep.cpp`. 
* To change the data structures that are tested in these experiments you can comment out the relevant lines in `benchmark/dynamic_trees/benchmark_updates_diameter_sweep.cpp` and `benchmark/dynamic_trees/benchmark_queries_diameter_sweep.cpp`. To save the most time we recommend not running the RC-Tree and Topology Tree baselines as they are significantly slower on some inputs.


### Parallel Update Speed Experiments (Figure 7)

This experiment uses all of the cores available by default. In our experiments we used 192 hyperthreads (see below for how to customize the number of threads). Run the shell script:

```bash
scripts/par_update_experiments.sh
```

* Results written to: `results/parallel_update_speed_10000000_1000000.csv` and `results/parallel_update_speed_graph.csv`.
* Results plotted in: `results/par_update_speed.pdf`.

Customize This Experiment:
* To customize the number of threads change the line that sets the environment variable `PARLAY_NUM_THREADS` in `scripts/par_update_experiments.sh`.
* These experiments take the average time over $3$ trials for each case. You can change this value by modifying the $4$-th value in each tuple entry in the `test_cases` lists in `benchmark/parallel/parallel_benchmark_runner.cpp` and `benchmark/graph_algorithms/benchmark_graph_update_speed_parallel.cpp`. To save time we recommend setting this to $1$ for all cases.
* You can change the values of `N` and `K` in `scripts/par_update_experiments.sh`.
  * `N` is the size of the synthetic trees for the parallel update speed experiment.
  * `K` is the batch size for both the synthetic and real-world experiments.
* To change the data structures that are tested in these experiments you can comment out the relevant lines in `benchmark/parallel/parallel_benchmark_runner.cpp` and `benchmark/graph_algorithms/benchmark_graph_update_speed_parallel.cpp`. To save the most time we recommend not running the RC-Tree and Topology Tree baselines as they are significantly slower on some inputs.


### Input Size Scaling Experiment (Figure 8)

This experiment uses all of the cores available by default. In our experiments we used 192 hyperthreads (see below for how to customize the number of threads). Run the shell script:

```bash
scripts/size_sweep_experiments.sh
```

* Results written to: `results/n_sweep_parallel_update.csv`.
* Results plotted in: `results/n_sweep.pdf`.

Customize This Experiment:
* To customize the number of threads change the line that sets the environment variable `PARLAY_NUM_THREADS` in `scripts/par_update_experiments.sh`.
* These experiments take the average time over $1$ trials for each case. You can change this value by modifying the $4$-th value in each tuple entry in the `test_cases` lists in `benchmark/parallel/parallel_benchmark_runner.cpp` and `benchmark/graph_algorithms/benchmark_graph_update_speed_parallel.cpp`. To save time we recommend leaving this at $1$ for all cases.
* To change the values of $n$ used in the sweep, modify the `n_list` variable in `benchmark/parallel/parallel_benchmark_n_sweep.cpp`.
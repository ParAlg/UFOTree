#pragma once
#include "../dynamic_trees/benchmark.h"
#include "parallel_cluster_benchmark.h"
#include "parallel_ufo_cluster_unordered_set.h"
#include "parallel_ufo_cluster_abseil_set.h"
#include "parallel_ufo_cluster_PAM.h"
#include "util.h"
using ll = long long;

vector<ll> batch_sizes({1,10,100,1000,10000});
vector<ll> granularities({1,10,100,1000,10000});

int main(){
  /*std::cout << "[Benchmarking Single Insert Delete Continuous:] \n\n";

  std::cout << "Unordered Set: " << "\n";
  benchmark_single_insert_delete_continuous<ParallelUFOClusterUSet<int>>();
  std::cout << "\n";
  std::cout << "Abseil Hash Set: " << "\n";
  benchmark_single_insert_delete_continuous<ParallelUFOClusterASet<int>>();
  std::cout << "\n";
  std::cout << "PAM: " << "\n";
  benchmark_single_insert_delete_continuous<ParallelUFOClusterPAM<int>>();
  std::cout << "\n";
  std::cout << "Other Parallel Hash Map: " << "\n";*/
  for(auto b_size : batch_sizes){
    batch_size = b_size;
    for(auto g : granularities){
      if(g > b_size) break;
      granularity = g; 
      std::cout << "[Benchmarking Single Insert Delete Randomized : Batch Size = " << batch_size << ", Granularity = " << granularity << "]\n\n";

      std::cout << "Unordered Set: " << "\n";
      benchmark_single_insert_delete_randomized<ParallelUFOClusterUSet<int>>();
      std::cout << "\n";
      std::cout << "Abseil Hash Set: " << "\n";
      benchmark_single_insert_delete_randomized<ParallelUFOClusterASet<int>>();
      std::cout << "PAM: " << "\n";
      benchmark_single_insert_delete_randomized<ParallelUFOClusterPAM<int>>();
      //std::cout << "Other Parallel Hash Map: " << "\n";

      std::cout << "[Benchmarking Batched Insert Delete: Batch Size = " << batch_size << ", Granularity = " << granularity << "]\n\n";
      std::cout << "Unordered Set: " << "\n";
      benchmark_batch_updates<ParallelUFOClusterUSet<int>>();
      std::cout << "\n";
      std::cout << "Abseil Hash Set: " << "\n";
      benchmark_batch_updates<ParallelUFOClusterASet<int>>();
      std::cout << "PAM: " << "\n";
      benchmark_batch_updates<ParallelUFOClusterPAM<int>>();
      std::cout << "\n";
      //std::cout << "Other Parallel Hash Map: " << "\n";
    }
  }
}


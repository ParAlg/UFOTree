#pragma once
#include "../dynamic_trees/benchmark.h"
#include "parallel_cluster_benchmark.h"
#include "parallel_ufo_cluster_unordered_set.h"
#include "parallel_ufo_cluster_abseil_set.h"
#include "parallel_ufo_cluster_PAM.h"
#include "util.h"
using ll = long long;

vector<ll> batch_sizes({10000000});
vector<ll> granularities({100});

std::tuple<std::string, std::function<std::vector<Update>(vertex_t, long)>, bool, int> test_cases[] = {
    {"Linked List", dynamic_tree_benchmark::linked_list_benchmark, false, 1},
    {"Binary Tree", dynamic_tree_benchmark::binary_tree_benchmark, false, 1},
    {"64-ary Tree", dynamic_tree_benchmark::k_ary_tree_benchmark, true, 1},
    {"Star", dynamic_tree_benchmark::star_benchmark, true, 1},
    {"Random Degree 3", dynamic_tree_benchmark::random_degree3_benchmark, false, 5},
    {"Random Unbounded Degree", dynamic_tree_benchmark::random_unbounded_benchmark, true, 5},
    {"Preferential Attachment", dynamic_tree_benchmark::preferential_attachment_benchmark, true, 5}
  };


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
      /*granularity = g; 
      std::cout << "[Benchmarking Single Insert Delete Randomized : Batch Size = " << batch_size << ", Granularity = " << granularity << "]\n\n";

      std::cout << "Unordered Set: " << "\n";
      benchmark_single_insert_delete_randomized<ParallelUFOClusterUSet<int>>();
      std::cout << "\n";
      std::cout << "Abseil Hash Set: " << "\n";
      benchmark_single_insert_delete_randomized<ParallelUFOClusterASet<int>>();
      std::cout << "PAM: " << "\n";
      //benchmark_single_insert_delete_randomized<ParallelUFOClusterPAM<int>>();
      //std::cout << "Other Parallel Hash Map: " << "\n";

      std::cout << "[Benchmarking Batched Insert Delete: Batch Size = " << batch_size << ", Granularity = " << granularity << "]\n\n";
      std::cout << "Unordered Set: " << "\n";
      //benchmark_batch_updates<ParallelUFOClusterUSet<int>>();
      std::cout << "\n";
      std::cout << "Abseil Hash Set: " << "\n";
      //benchmark_batch_updates<ParallelUFOClusterASet<int>>();
      std::cout << "PAM: " << "\n";
      benchmark_batch_updates<ParallelUFOClusterPAM<int>>();
      std::cout << "\n";
      //std::cout << "Other Parallel Hash Map: " << "\n";*/

      std::cout << "[Benchmarking Real Trees with Batch Size = " << batch_size << ", Granularity = " << granularity << "]\n\n";
      std::cout << "PAM: " << "\n";
      for(auto test_case: test_cases){
        std::cout << "Tree Type: " << std::get<0>(test_case) << "\n";
        benchmark_real_trees<ParallelUFOClusterPAM<int>>(batch_size, std::get<1>(test_case));
        std::cout << "\n";
      }
  }
}
}

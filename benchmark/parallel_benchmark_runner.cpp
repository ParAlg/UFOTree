#include "benchmark.h"
#include "parallel_benchmark.h"
#include "util.h"
#include "parallel_ufo_tree.h"
#include "parallel_topology_tree.h"
// #include "../baselines/dynamic_trees/parallel_euler_tour_tree/include/euler_tour_tree.hpp"
#include <fstream>


int main(int argc, char** argv) {
  // List of values of n to loop through and run all test cases
  std::vector<vertex_t> n_list = {1000};
  if (argc < 2) {
    std::cout << "Using default hard-coded list for values of n." << std::endl;
  } else {
    std::cout << "Using command line arguments for values of n." << std::endl;
    n_list.clear();
    for (int i = 1; i < argc; ++i) n_list.push_back(std::atoi(argv[i]));
  }
  vertex_t k = 100;
  /* Each test case has a name for output, the update generator function, and
  a bool indicating if ternarization may be necessary for this input */
  std::tuple<std::string, std::function<std::vector<Update>(vertex_t)>, bool, int> test_cases[] = {
    {"Linked List", dynamic_tree_benchmark::linked_list_benchmark, false, 1},
    {"Binary Tree", dynamic_tree_benchmark::binary_tree_benchmark, false, 1},
    {"64-ary Tree", dynamic_tree_benchmark::k_ary_tree_benchmark, true, 1},
    {"Star", dynamic_tree_benchmark::star_benchmark, true, 1},
    {"Random Degree 3", dynamic_tree_benchmark::random_degree3_benchmark, false, 5},
    {"Random Unbounded Degree", dynamic_tree_benchmark::random_unbounded_benchmark, true, 5},
    {"Preferential Attachment", dynamic_tree_benchmark::preferential_attachment_benchmark, true, 5}
  };

  for (vertex_t n : n_list) {
    std::string filename = "../results/parellel_update_speed_" + std::to_string(n) + "_" + std::to_string(k) + ".csv";
    std::ofstream output_csv;
    output_csv.open(filename);
    output_csv << "Test Case,Euler Tour Tree,UFO Tree,Topology Tree,\n";

    for (auto test_case : test_cases) {
      std::string test_case_name = std::get<0>(test_case);
      auto update_generator = std::get<1>(test_case);
      bool ternarize = std::get<2>(test_case);
      int num_trials = std::get<3>(test_case);
      std::vector<std::vector<Update>> update_sequences;
      for (int i = 0; i < num_trials; i++) {
        std::vector<Update> updates = update_generator(n);
        update_sequences.push_back(updates);
      }
      double time;
      std::cout << "[ RUNNING " << test_case_name << " PARALLEL UPDATE SPEED BENCHMARK WITH n=" << n << " ]" << std::endl;
      output_csv << test_case_name << ",";

      // Euler Tour Tree
      // time = parallel_dynamic_tree_benchmark::get_update_speed<ParallelEulerTourTree>(n, k, update_sequences);
      // std::cout << "EulerTourTree : " << time << std::endl;
      // output_csv << time << ",";
      // UFO Tree
      time = parallel_dynamic_tree_benchmark::get_update_speed<ParallelUFOTree<int>>(n, k, update_sequences);
      std::cout << "UFOTree       : " << time << std::endl;
      output_csv << time << ",";
      // Topology Tree
      if (!ternarize) {
        time = parallel_dynamic_tree_benchmark::get_update_speed<ParallelTopologyTree<int>>(n, k, update_sequences);
        std::cout << "TopologyTree  : " << time << std::endl;
        output_csv << time << ",";
      }

      std::cout << std::endl;
      output_csv << "\n";
    }

    output_csv.close();
  }
}

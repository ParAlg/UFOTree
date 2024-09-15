#include "benchmark.h"
#include "../include/util.h"
#include "../include/ufo_tree.h"
#include "../include/topology_tree.h"
#include "../include/rc_tree.h"
#include "../baselines/dynamic_trees/euler_tour_tree/include/skip_list_ett.hpp"
#include <fstream>


using namespace skip_list_ett;

int main(int argc, char** argv) {
  // List of values of n to loop through and run all test cases
  vertex_t n_list[] = {1000};
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
    std::string filename = "../results/update_speed_" + std::to_string(n) + ".csv";
    std::ofstream output_csv;
    output_csv.open(filename);
    output_csv << "Test Case,RC Tree,Topology Tree,UFO Tree,Euler Tour Tree,\n";

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
      std::cout << "[ RUNNING " << test_case_name << " BENCHMARK WITH n=" << n << " ]" << std::endl;
      output_csv << test_case_name << ",";

      // RC Tree
      if (!ternarize) time = dynamic_tree_benchmark::get_update_speed<RCTree<int>>(n, update_sequences);
      else time = 0;
      std::cout << "RCTree        : " << time << std::endl;
      output_csv << time << ",";
      // Topology Tree
      if (!ternarize) time = dynamic_tree_benchmark::get_update_speed<TopologyTree<int>>(n, update_sequences);
      else time = 0;
      std::cout << "TopologyTree  : " << time << std::endl;
      output_csv << time << ",";
      // UFO Tree
      time = dynamic_tree_benchmark::get_update_speed<UFOTree<int>>(n, update_sequences);
      std::cout << "UFOTree       : " << time << std::endl;
      output_csv << time << ",";
      // Euler Tour Tree
      time = dynamic_tree_benchmark::get_update_speed<EulerTourTree>(n, update_sequences);
      std::cout << "EulerTourTree : " << time << std::endl;
      output_csv << time << ",";

      std::cout << std::endl;
      output_csv << "\n";
    }

    output_csv.close();
  }
}

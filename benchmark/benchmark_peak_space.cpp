#include "benchmark.h"
#include "util.h"
#include "ternarized_tree.h"
#include "ufo_tree.h"
#include "topology_tree.h"
#include "rc_tree.h"
#include "../baselines/dynamic_trees/euler_tour_tree/include/skip_list_ett.hpp"
#include <fstream>


using namespace skip_list_ett;

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
  srand(time(NULL));
  /* Each test case has a name for output, the update generator function, and
  a bool indicating if ternarization may be necessary for this input */
  std::tuple<std::string, std::function<std::vector<Update>(vertex_t, long)>, bool, int> test_cases[] = {
    {"Linked List", dynamic_tree_benchmark::linked_list_benchmark, false, 1},
    {"Binary Tree", dynamic_tree_benchmark::binary_tree_benchmark, false, 1},
    {"64-ary Tree", dynamic_tree_benchmark::k_ary_tree_benchmark, true, 1},
    {"Star", dynamic_tree_benchmark::star_benchmark, true, 1},
    {"Dandelion", dynamic_tree_benchmark::dandelion_benchmark, true, 1},
    {"Random Degree 3", dynamic_tree_benchmark::random_degree3_benchmark, false, 1},
    {"Random Unbounded Degree", dynamic_tree_benchmark::random_unbounded_benchmark, true, 1},
    {"Preferential Attachment", dynamic_tree_benchmark::preferential_attachment_benchmark, true, 1}
  };

  for (vertex_t n : n_list) {
    std::string filename = "../results/peak_space_" + std::to_string(n) + ".csv";
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
        std::vector<Update> updates = update_generator(n, rand());
        update_sequences.push_back(updates);
      }
      double space;
      std::cout << "[ RUNNING " << test_case_name << " PEAK SPACE BENCHMARK WITH n=" << n << " ]" << std::endl;
      output_csv << test_case_name << ",";

      // UFO Tree
      space = dynamic_tree_benchmark::get_peak_space<UFOTree<int, empty_t>>(n, update_sequences);
      std::cout << "UFOTree       : " << space << std::endl;
      output_csv << space << ",";
      // Euler Tour Tree
      space = dynamic_tree_benchmark::get_peak_space<EulerTourTree>(n, update_sequences);
      std::cout << "EulerTourTree : " << space << std::endl;
      output_csv << space << ",";
      // RC Tree
      if (!ternarize) space = dynamic_tree_benchmark::get_peak_space<RCTree<int>>(n, update_sequences);
      else space = dynamic_tree_benchmark::get_peak_space<TernarizedTree<RCTree<int>, int>>(n, update_sequences);
      std::cout << "RCTree        : " << space << std::endl;
      output_csv << space << ",";
      // Topology Tree
      if (!ternarize) space = dynamic_tree_benchmark::get_peak_space<TopologyTree<int, empty_t>>(n, update_sequences);
      else space = dynamic_tree_benchmark::get_peak_space<TernarizedTree<TopologyTree<int, empty_t>, empty_t>>(n, update_sequences);
      std::cout << "TopologyTree  : " << space << std::endl;
      output_csv << space << ",";

      std::cout << std::endl;
      output_csv << "\n";
    }

    output_csv.close();
  }
}

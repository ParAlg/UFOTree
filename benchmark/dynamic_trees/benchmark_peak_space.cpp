#include "benchmark.h"
#include "util.h"
#include "ternarized_tree.h"
#include "ufo_tree.h"
#include "topology_tree.h"
#include "rc_tree.h"
#include "spaa_rc_tree.h"
#include "spaa_rc_tree_ternarized.h"
#include "parett/dynamic_trees/link_cut_tree/link_cut_tree.hpp"
#include "parett/dynamic_trees/euler_tour_tree/skip_list_ett.hpp"
#include "parett/dynamic_trees/euler_tour_tree/splay_tree_ett.hpp"
#include "parett/dynamic_trees/euler_tour_tree/treap_ett.hpp"
#include "top_tree.h"
#include <fstream>

using namespace ufo;


int main(int argc, char** argv) {
  // List of values of n to loop through and run all test cases
  std::vector<vertex_t> n_list = {10000000};
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
    //output_csv << "Test Case,UFO Tree,Euler Tour Tree,RC Tree,Topology Tree,\n";
    output_csv << "Test Case, UFO Tree, ETT (Skip List), Splay Top Tree, ETT (Splay Tree), ETT (Treap), Link-Cut Tree";
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
      space = dynamic_tree_benchmark::get_peak_space<skip_list_ett::EulerTourTree>(n, update_sequences);
      std::cout << "EulerTourTree : " << space << std::endl;
      output_csv << space << ",";
      /*// RC Tree
      if (!ternarize) space = dynamic_tree_benchmark::get_peak_space<RCTree<int>>(n, update_sequences);
      else space = dynamic_tree_benchmark::get_peak_space<TernarizedTree<RCTree<int>, int>>(n, update_sequences);
      std::cout << "RCTree        : " << space << std::endl;
      output_csv << space << ",";

      // Topology Tree
      if (!ternarize) space = dynamic_tree_benchmark::get_peak_space<TopologyTree<int, empty_t>>(n, update_sequences);
      else space = dynamic_tree_benchmark::get_peak_space<TernarizedTree<TopologyTree<int, empty_t>, empty_t>>(n, update_sequences);
      std::cout << "TopologyTree  : " << space << std::endl;
      output_csv << space << ",";

      */
      
      // Splay Top Tree
      space = dynamic_tree_benchmark::get_peak_space<TopTree<int>>(n, update_sequences);
      std::cout << "Splay Top Tree: " << space << std::endl;
      output_csv << space << ",";

      // Splay Euler Tour Tree
      space = dynamic_tree_benchmark::get_peak_space<splay_tree_ett::EulerTourTree>(n, update_sequences);
      std::cout << "EulerTourTree (Splay Tree): " << space << std::endl;
      output_csv << space << ",";


      // Treap Euler Tour Tree
      space = dynamic_tree_benchmark::get_peak_space<treap::EulerTourTree<int>>(n, update_sequences);
      std::cout << "EulerTourTree (Treap): " << space << std::endl;
      output_csv << space << ",";


      // Link-Cut Tree  
      space = dynamic_tree_benchmark::get_peak_space<link_cut_tree::LinkCutTreeInt>(n, update_sequences);
      std::cout << "Link Cut Tree: " << space << std::endl;
      output_csv << space;



      // CMU RC Tree
      /*if (!ternarize) space = dynamic_tree_benchmark::get_peak_space<ParallelRCTree<int>>(n, update_sequences);
      else space = dynamic_tree_benchmark::get_peak_space<ParallelRCTreeTernarized<int>>(n, update_sequences);
      std::cout << "RCTree (CMU Version)     : " << space << std::endl;
      if(test_case_name != "Preferential Attachment") output_csv<< ",";*/


      std::cout << std::endl;
      output_csv << "\n";
    }

    output_csv.close();
  }
}

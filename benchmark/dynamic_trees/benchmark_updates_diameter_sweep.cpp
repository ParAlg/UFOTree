#include "benchmark.h"
#include "util.h"
#include "ternarized_tree.h"
#include "ufo_tree.h"
#include "topology_tree.h"
#include "rc_tree.h"
#include "top_tree.h"
#include "parett/dynamic_trees/link_cut_tree/link_cut_tree.hpp"
#include "parett/dynamic_trees/euler_tour_tree/skip_list_ett.hpp"
#include "parett/dynamic_trees/euler_tour_tree/splay_tree_ett.hpp"
#include "parett/dynamic_trees/euler_tour_tree/treap_ett.hpp"
#include <fstream>

using namespace ufo;


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
    {"0.00", dynamic_tree_benchmark::zipf_tree_benchmark<0.00>, true, 3},
    {"0.50", dynamic_tree_benchmark::zipf_tree_benchmark<0.50>, true, 3},
    {"1.01", dynamic_tree_benchmark::zipf_tree_benchmark<1.01>, true, 3},
    {"1.50", dynamic_tree_benchmark::zipf_tree_benchmark<1.50>, true, 3},
    {"2.00", dynamic_tree_benchmark::zipf_tree_benchmark<2.00>, true, 3},
  };

  for (vertex_t n : n_list) {
    // Update speed while supporting queries benchmark
    // std::string filename = "../results/diameter_sweep_update_speed_" + std::to_string(n) + ".csv";
    // std::ofstream output_csv;
    // output_csv.open(filename);
    // output_csv << "Test Case,"
    //             << "Link-Cut Tree,"
    //             << "UFO Tree,"
    //             << "Splay Top Tree,"
    //             << "ETT (Treap),"
    //             << "ETT (Splay Tree),"
    //             << "ETT (Skip List),"
    //             << "Topology Tree,"
    //             << "Rake-Compress Tree";

    // for (auto test_case : test_cases) {
    //   std::string test_case_name = std::get<0>(test_case);
    //   auto update_generator = std::get<1>(test_case);
    //   bool ternarize = std::get<2>(test_case);
    //   int num_trials = std::get<3>(test_case);
    //   std::vector<std::vector<Update>> update_sequences;
    //   for (int i = 0; i < num_trials; i++) {
    //     std::vector<Update> updates = update_generator(n, rand());
    //     // std::vector<Update> updates = dynamic_tree_benchmark::dense_tree_update_generator(update_generator(n, rand()));
    //     update_sequences.push_back(updates);
    //   }
    //   double time;
    //   std::cout << "[ RUNNING " << test_case_name << " INT UPDATE SPEED BENCHMARK WITH n=" << n << " trials=" << num_trials << " ]" << std::endl;
    //   output_csv << "\n" << test_case_name;

    //   // Link Cut Tree
    //   time = dynamic_tree_benchmark::get_update_speed<link_cut_tree::LinkCutTreeInt>(n, update_sequences);
    //   std::cout << "LinkCutTree   : " << time << std::endl;
    //   output_csv << "," << time;
    //   // UFO Tree
    //   time = dynamic_tree_benchmark::get_update_speed<UFOTree<int, int>>(n, update_sequences);
    //   std::cout << "UFOTree       : " << time << std::endl;
    //   output_csv << "," << time;
    //   // Top Tree
    //   time = dynamic_tree_benchmark::get_update_speed<TopTree<int>>(n, update_sequences);
    //   std::cout << "SplayTopTree  : " << time << std::endl;
    //   output_csv << "," << time;
    //   // Euler Tour Tree (Treap)
    //   time = dynamic_tree_benchmark::get_update_speed<treap::EulerTourTree<int>>(n, update_sequences);
    //   std::cout << "TreapETT      : " << time << std::endl;
    //   output_csv << "," << time;
    //   // Euler Tour Tree (Splay Tree)
    //   time = dynamic_tree_benchmark::get_update_speed<splay_tree_ett::EulerTourTree>(n, update_sequences);
    //   std::cout << "SplayTreeETT  : " << time << std::endl;
    //   output_csv << "," << time;
    //   // Euler Tour Tree (Skip List)
    //   time = dynamic_tree_benchmark::get_update_speed<skip_list_ett::EulerTourTree>(n, update_sequences);
    //   std::cout << "SkipListETT   : " << time << std::endl;
    //   output_csv << "," << time;
    //   // Topology Tree
    //   if (!ternarize) time = dynamic_tree_benchmark::get_update_speed<TopologyTree<int, int>>(n, update_sequences);
    //   else time = dynamic_tree_benchmark::get_update_speed<TernarizedTree<TopologyTree<int, int>, int>>(n, update_sequences);
    //   std::cout << "TopologyTree  : " << time << std::endl;
    //   output_csv << "," << time;
    //   // RC Tree
    //   if (!ternarize) time = dynamic_tree_benchmark::get_update_speed<RCTree<int>>(n, update_sequences);
    //   else time = dynamic_tree_benchmark::get_update_speed<TernarizedTree<RCTree<int>, int>>(n, update_sequences);
    //   std::cout << "RCTree        : " << time << std::endl;
    //   output_csv << "," << time;

    //   std::cout << std::endl;
    // }
    // output_csv.close();

    // Empty dynamic tree update speed benchmark
    std::string filename_empty = "../results/diameter_sweep_update_" + std::to_string(n) + ".csv";
    std::ofstream output_csv_empty;
    output_csv_empty.open(filename_empty);
    output_csv_empty << "Alpha,"
                << "Link-Cut Tree,"
                << "UFO Tree,"
                << "Splay Top Tree,"
                << "ETT (Treap),"
                << "ETT (Splay Tree),"
                << "ETT (Skip List),"
                << "Topology Tree,"
                << "Rake-Compress Tree";

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
      double time;
      std::cout << "[ RUNNING " << test_case_name << " EMPTY UPDATE SPEED BENCHMARK WITH n=" << n << " trials=" << num_trials << " ]" << std::endl;
      output_csv_empty << "\n" << test_case_name;

      // Link Cut Tree
      time = dynamic_tree_benchmark::get_update_speed<link_cut_tree::LinkCutTree>(n, update_sequences);
      std::cout << "LinkCutTree   : " << time << std::endl;
      output_csv_empty << "," << time;
      // UFO Tree
      time = dynamic_tree_benchmark::get_update_speed<UFOTree<empty_t, empty_t>>(n, update_sequences);
      std::cout << "UFOTree       : " << time << std::endl;
      output_csv_empty << "," << time;
      // Top Tree
      time = dynamic_tree_benchmark::get_update_speed<TopTree<empty_t>>(n, update_sequences);
      std::cout << "SplayTopTree  : " << time << std::endl;
      output_csv_empty << "," << time;
      // Euler Tour Tree (Treap)
      time = dynamic_tree_benchmark::get_update_speed<treap::EulerTourTree<empty_t>>(n, update_sequences);
      std::cout << "TreapETT      : " << time << std::endl;
      output_csv_empty << "," << time;
      // Euler Tour Tree (Splay Tree)
      time = dynamic_tree_benchmark::get_update_speed<splay_tree_ett::EulerTourTree>(n, update_sequences);
      std::cout << "SplayTreeETT  : " << time << std::endl;
      output_csv_empty << "," << time;
      // Euler Tour Tree (Skip List)
      time = dynamic_tree_benchmark::get_update_speed<skip_list_ett::EulerTourTree>(n, update_sequences);
      std::cout << "SkipListETT   : " << time << std::endl;
      output_csv_empty << "," << time;
      // Topology Tree
      if (!ternarize) time = dynamic_tree_benchmark::get_update_speed<TopologyTree<empty_t, empty_t>>(n, update_sequences);
      else time = dynamic_tree_benchmark::get_update_speed<TernarizedTree<TopologyTree<empty_t, empty_t>, empty_t>>(n, update_sequences);
      std::cout << "TopologyTree  : " << time << std::endl;
      output_csv_empty << "," << time;
      // RC Tree
      if (!ternarize) time = dynamic_tree_benchmark::get_update_speed<RCTree<int>>(n, update_sequences);
      else time = dynamic_tree_benchmark::get_update_speed<TernarizedTree<RCTree<int>, int>>(n, update_sequences);
      std::cout << "RCTree        : " << time << std::endl;
      output_csv_empty << "," << time;

      std::cout << std::endl;
    }
    output_csv_empty.close();
  }
}

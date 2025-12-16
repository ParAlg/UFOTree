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
    {"Linked List", dynamic_tree_benchmark::linked_list_benchmark, false, 1},
    {"Binary Tree", dynamic_tree_benchmark::binary_tree_benchmark, false, 1},
    {"64-ary Tree", dynamic_tree_benchmark::k_ary_tree_benchmark, true, 1},
    {"Star", dynamic_tree_benchmark::star_benchmark, true, 1},
    {"Dandelion", dynamic_tree_benchmark::dandelion_benchmark, true, 1},
    {"Random Degree 3", dynamic_tree_benchmark::random_degree3_benchmark, false, 1},
    {"Random Unbounded Degree", dynamic_tree_benchmark::random_unbounded_benchmark, true, 1},
    {"Preferential Attachment", dynamic_tree_benchmark::preferential_attachment_benchmark, true, 1},
  };

  for (vertex_t n : n_list) {
    std::string filename = "../results/query_speed_" + std::to_string(n) + ".csv";
    std::ofstream output_csv;
    output_csv.open(filename);

    // Connectivity query speed benchmark
    output_csv << "Test Case,UFO Tree,Link Cut Tree,Splay Top Tree,ETT (Skip List),ETT (Splay Tree),ETT (Treap),Topology Tree,RC Tree";
    for (auto test_case : test_cases) {
      std::string test_case_name = std::get<0>(test_case);
      auto update_generator = std::get<1>(test_case);
      bool ternarize = std::get<2>(test_case);
      int num_trials = std::get<3>(test_case);
      std::vector<std::vector<Update>> update_sequences;
      std::vector<std::vector<Query>> query_sequences;
      for (int i = 0; i < num_trials; i++) {
        std::vector<Update> updates = update_generator(n, rand());
        update_sequences.push_back(updates);
        std::vector<Query> queries = dynamic_tree_benchmark::random_query_generator(n, n);
        query_sequences.push_back(queries);
      }
      double time;
      std::cout << "[ RUNNING " << test_case_name << " CONNECTIVITY QUERY SPEED BENCHMARK WITH n=" << n << " ]" << std::endl;
      output_csv << "\n" << test_case_name << "_conn";

      // UFO Tree
      time = dynamic_tree_benchmark::get_conn_query_speed<UFOTree<empty_t, empty_t>>(n, update_sequences, query_sequences);
      std::cout << "UFOTree       : " << time << std::endl;
      output_csv << "," << time;
      // Link Cut Tree
      time = dynamic_tree_benchmark::get_conn_query_speed<link_cut_tree::LinkCutTree>(n, update_sequences, query_sequences);
      std::cout << "LinkCutTree   : " << time << std::endl;
      output_csv << "," << time;
      // Top Tree
      time = dynamic_tree_benchmark::get_conn_query_speed<TopTree<empty_t>>(n, update_sequences, query_sequences);
      std::cout << "SplayTopTree  : " << time << std::endl;
      output_csv << "," << time;
      // Euler Tour Tree (Skip List)
      time = dynamic_tree_benchmark::get_conn_query_speed<skip_list_ett::EulerTourTree>(n, update_sequences, query_sequences);
      std::cout << "SkipListETT   : " << time << std::endl;
      output_csv << "," << time;
      // Euler Tour Tree (Splay Tree)
      time = dynamic_tree_benchmark::get_conn_query_speed<splay_tree_ett::EulerTourTree>(n, update_sequences, query_sequences);
      std::cout << "SplayTreeETT  : " << time << std::endl;
      output_csv << "," << time;
      // Euler Tour Tree (Treap)
      time = dynamic_tree_benchmark::get_conn_query_speed<treap::EulerTourTree<empty_t>>(n, update_sequences, query_sequences);
      std::cout << "TreapETT      : " << time << std::endl;
      output_csv << "," << time;
      // Topology Tree
      if (!ternarize) time = dynamic_tree_benchmark::get_conn_query_speed<TopologyTree<int, empty_t>>(n, update_sequences, query_sequences);
      else time = dynamic_tree_benchmark::get_conn_query_speed<TernarizedTree<TopologyTree<int, empty_t>, empty_t>>(n, update_sequences, query_sequences);
      std::cout << "TopologyTree  : " << time << std::endl;
      output_csv << "," << time;
      // RC Tree
      if (!ternarize) time = dynamic_tree_benchmark::get_conn_query_speed<RCTree<int>>(n, update_sequences, query_sequences);
      else time = dynamic_tree_benchmark::get_conn_query_speed<TernarizedTree<RCTree<int>, int>>(n, update_sequences, query_sequences);
      std::cout << "RCTree        : " << time << std::endl;
      output_csv << "," << time;

      std::cout << std::endl;
    }

    // Path query speed benchmark
    output_csv << "\n";
    output_csv << "Test Case,UFO Tree,Link Cut Tree,Splay Top Tree,Topology Tree,RC Tree";
    for (auto test_case : test_cases) {
      std::string test_case_name = std::get<0>(test_case);
      auto update_generator = std::get<1>(test_case);
      bool ternarize = std::get<2>(test_case);
      int num_trials = std::get<3>(test_case);
      std::vector<std::vector<Update>> update_sequences;
      std::vector<std::vector<Query>> query_sequences;
      for (int i = 0; i < num_trials; i++) {
        std::vector<Update> updates = update_generator(n, rand());
        update_sequences.push_back(updates);
        std::vector<Query> queries = dynamic_tree_benchmark::random_query_generator(n, n);
        query_sequences.push_back(queries);
      }
      double time;
      std::cout << "[ RUNNING " << test_case_name << " PATH QUERY SPEED BENCHMARK WITH n=" << n << " ]" << std::endl;
      output_csv << "\n" << test_case_name << "_path";

      // UFO Tree
      time = dynamic_tree_benchmark::get_path_query_speed<UFOTree<int, int>>(n, update_sequences, query_sequences);
      std::cout << "UFOTree       : " << time << std::endl;
      output_csv << "," << time;
      // Link Cut Tree
      time = dynamic_tree_benchmark::get_path_query_speed<link_cut_tree::LinkCutTreeInt>(n, update_sequences, query_sequences);
      std::cout << "LinkCutTree   : " << time << std::endl;
      output_csv << "," << time;
      //Top Tree
      time = dynamic_tree_benchmark::get_path_query_speed<TopTree<int>>(n, update_sequences, query_sequences);
      std::cout << "SplayTopTree  : " << time << std::endl;
      output_csv << "," << time;
      // Topology Tree
      if (!ternarize) time = dynamic_tree_benchmark::get_path_query_speed<TopologyTree<int, int>>(n, update_sequences, query_sequences);
      else time = dynamic_tree_benchmark::get_path_query_speed<TernarizedTree<TopologyTree<int, int>, int>>(n, update_sequences, query_sequences);
      std::cout << "TopologyTree  : " << time << std::endl;
      output_csv << "," << time;
      // RC Tree
      if (!ternarize) time = dynamic_tree_benchmark::get_path_query_speed<RCTree<int>>(n, update_sequences, query_sequences);
      else time = dynamic_tree_benchmark::get_path_query_speed<TernarizedTree<RCTree<int>, int>>(n, update_sequences, query_sequences);
      std::cout << "RCTree        : " << time << std::endl;
      output_csv << "," << time;

      std::cout << std::endl;
    }

    output_csv.close();
  }
}

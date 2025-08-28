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

using namespace dgbs;


int main(int argc, char** argv) {
  // List of values of n to loop through and run all test cases
  vertex_t n = 100000;
  vertex_t q = 1000000;
  if (argc < 2) {
      std::cout << "Using default hard-coded list for values of n and q." << std::endl;
  } else if (argc == 3) {
      std::cout << "Using command line arguments for values of n and q." << std::endl;
      n = std::stoi(argv[1]);
      q = std::stoi(argv[2]);
  } else {
      std::cout << "Usage: ./benchmark [n] [q]" << std::endl;
  }
  long seed = 0;
  srand(0);
  /* Each test case has a name for output, the update generator function, and
  a bool indicating if ternarization may be necessary for this input */
  std::tuple<std::string, std::function<std::vector<Update>(vertex_t, long)>, bool, int> test_cases[] = {
    {"0.00", dynamic_tree_benchmark::zipf_tree_benchmark<0.00>, true, 3},
    // {"0.25", dynamic_tree_benchmark::zipf_tree_benchmark<0.25>, true, 3},
    {"0.50", dynamic_tree_benchmark::zipf_tree_benchmark<0.50>, true, 3},
    // {"0.75", dynamic_tree_benchmark::zipf_tree_benchmark<0.75>, true, 3},
    {"1.01", dynamic_tree_benchmark::zipf_tree_benchmark<1.01>, true, 3},
    // {"1.25", dynamic_tree_benchmark::zipf_tree_benchmark<1.25>, true, 3},
    {"1.50", dynamic_tree_benchmark::zipf_tree_benchmark<1.50>, true, 3},
    // {"1.75", dynamic_tree_benchmark::zipf_tree_benchmark<1.75>, true, 3},
    {"2.00", dynamic_tree_benchmark::zipf_tree_benchmark<2.00>, true, 3},
    // {"4.00", dynamic_tree_benchmark::zipf_tree_benchmark<4.00>, true, 3},
  };

  // Path query diameter sweep benchmark
  std::string filename_path = "../results/diameter_sweep_path_query_" + std::to_string(n) + ".csv";
  std::ofstream output_csv_path;
  output_csv_path.open(filename_path);
  output_csv_path << "Alpha,"
              << "Link-Cut Tree,"
              << "UFO Tree,"
              << "Splay Top Tree,"
              << "Topology Tree,"
              << "Rake-Compress Tree";

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
      std::vector<Query> queries = dynamic_tree_benchmark::random_query_generator(n, q);
      query_sequences.push_back(queries);
    }
    double time;
    std::cout << "[ RUNNING " << test_case_name << " PATH QUERY BENCHMARK WITH n=" << n << ", q=" << q << ", trials=" << num_trials << " ]" << std::endl;
    output_csv_path << "\n" << test_case_name;

    // Link Cut Tree
    time = dynamic_tree_benchmark::get_path_query_speed<link_cut_tree::LinkCutTreeInt>(n, update_sequences, query_sequences);
    std::cout << "LinkCutTree   : " << time << std::endl;
    output_csv_path << "," << time;
    // UFO Tree
    time = dynamic_tree_benchmark::get_path_query_speed<UFOTree<int, int>>(n, update_sequences, query_sequences);
    std::cout << "UFOTree       : " << time << std::endl;
    output_csv_path << "," << time;
    //Top Tree
    time = dynamic_tree_benchmark::get_path_query_speed<TopTree<int>>(n, update_sequences, query_sequences);
    std::cout << "SplayTopTree  : " << time << std::endl;
    output_csv_path << "," << time;
    // Topology Tree
    if (!ternarize) time = dynamic_tree_benchmark::get_path_query_speed<TopologyTree<int, int>>(n, update_sequences, query_sequences);
    else time = dynamic_tree_benchmark::get_path_query_speed<TernarizedTree<TopologyTree<int, int>, int>>(n, update_sequences, query_sequences);
    std::cout << "TopologyTree  : " << time << std::endl;
    output_csv_path << "," << time;
    // RC Tree
    if (!ternarize) time = dynamic_tree_benchmark::get_path_query_speed<RCTree<int>>(n, update_sequences, query_sequences);
    else time = dynamic_tree_benchmark::get_path_query_speed<TernarizedTree<RCTree<int>, int>>(n, update_sequences, query_sequences);
    std::cout << "RCTree        : " << time << std::endl;
    output_csv_path << "," << time;

    std::cout << std::endl;
  }
  output_csv_path.close();

  // Connectivity query diameter sweep benchmark
  std::string filename_conn = "../results/diameter_sweep_conn_query_" + std::to_string(n) + ".csv";
  std::ofstream output_csv_conn;
  output_csv_conn.open(filename_conn);
  output_csv_conn << "Alpha,"
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
    std::vector<std::vector<Query>> query_sequences;
    for (int i = 0; i < num_trials; i++) {
      std::vector<Update> updates = update_generator(n, rand());
      update_sequences.push_back(updates);
      std::vector<Query> queries = dynamic_tree_benchmark::random_query_generator(n, q);
      query_sequences.push_back(queries);
    }
    double time;
    std::cout << "[ RUNNING " << test_case_name << " CONN QUERY BENCHMARK WITH n=" << n << ", q=" << q << ", trials=" << num_trials << " ]" << std::endl;
    output_csv_conn << "\n" << test_case_name;

    // Link Cut Tree
    time = dynamic_tree_benchmark::get_conn_query_speed<link_cut_tree::LinkCutTree>(n, update_sequences, query_sequences);
    std::cout << "LinkCutTree   : " << time << std::endl;
    output_csv_conn << "," << time;
    // UFO Tree
    time = dynamic_tree_benchmark::get_conn_query_speed<UFOTree<empty_t, empty_t>>(n, update_sequences, query_sequences);
    std::cout << "UFOTree       : " << time << std::endl;
    output_csv_conn << "," << time;
    // Top Tree
    time = dynamic_tree_benchmark::get_conn_query_speed<TopTree<empty_t>>(n, update_sequences, query_sequences);
    std::cout << "SplayTopTree  : " << time << std::endl;
    output_csv_conn << "," << time;
    // Euler Tour Tree (Treap)
    time = dynamic_tree_benchmark::get_conn_query_speed<treap::EulerTourTree<empty_t>>(n, update_sequences, query_sequences);
    std::cout << "TreapETT      : " << time << std::endl;
    output_csv_conn << "," << time;
    // Euler Tour Tree (Splay Tree)
    time = dynamic_tree_benchmark::get_conn_query_speed<splay_tree_ett::EulerTourTree>(n, update_sequences, query_sequences);
    std::cout << "SplayTreeETT  : " << time << std::endl;
    output_csv_conn << "," << time;
    // Euler Tour Tree (Skip List)
    time = dynamic_tree_benchmark::get_conn_query_speed<skip_list_ett::EulerTourTree>(n, update_sequences, query_sequences);
    std::cout << "SkipListETT   : " << time << std::endl;
    output_csv_conn << "," << time;
    // Topology Tree
    if (!ternarize) time = dynamic_tree_benchmark::get_conn_query_speed<TopologyTree<empty_t, empty_t>>(n, update_sequences, query_sequences);
    else time = dynamic_tree_benchmark::get_conn_query_speed<TernarizedTree<TopologyTree<empty_t, empty_t>, empty_t>>(n, update_sequences, query_sequences);
    std::cout << "TopologyTree  : " << time << std::endl;
    output_csv_conn << "," << time;
    // RC Tree
    if (!ternarize) time = dynamic_tree_benchmark::get_conn_query_speed<RCTree<int>>(n, update_sequences, query_sequences);
    else time = dynamic_tree_benchmark::get_conn_query_speed<TernarizedTree<RCTree<int>, int>>(n, update_sequences, query_sequences);
    std::cout << "RCTree        : " << time << std::endl;
    output_csv_conn << "," << time;

    std::cout << std::endl;
  }
  output_csv_conn.close();

}

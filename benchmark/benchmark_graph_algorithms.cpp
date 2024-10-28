#include "benchmark.h"
#include "graph_benchmark.h"
#include "util.h"
#include "ternarized_tree.h"
#include "ufo_tree.h"
#include "topology_tree.h"
#include "rc_tree.h"
#include "../baselines/dynamic_trees/euler_tour_tree/include/skip_list_ett.hpp"
#include <fstream>

using namespace skip_list_ett;

int main(int argc, char** argv){

  srand(time(NULL));
  std::tuple<std::string, int> test_cases[] = {
    {"/ssd1/zhongqi/graphdata/sym/com-youtube_sym.bin", 1},
    {"/ssd1/zhongqi/graphdata/sym/com-orkut_sym.bin", 1},
    {"/ssd1/zhongqi/graphdata/sym/soc-LiveJournal1_sym.bin", 1},
    {"/ssd1/zhongqi/graphdata/sym/twitter_sym.bin", 1},
    {"/ssd1/zhongqi/graphdata/sym/friendster_sym.bin", 1},
  };

  std::string filename = "../results/update_speed_graph.csv";
  std::ofstream output_csv;
  output_csv.open(filename);
  output_csv << "Test Case,RC Tree,Topology Tree,UFO Tree,Euler Tour Tree,\n";


  // Incremental Minimum Spanning Tree
  for (auto test_case: test_cases) {
    std::string file_name = std::get<0>(test_case);
    int num_trials = std::get<1>(test_case);
    auto G = graph_utils::break_sym_graph_from_bin(file_name);
    auto n = G.size();
    parlay::sequence<std::pair<int, Edge>> updates;
    for (int i = 0; i < num_trials; i++) {
      auto non_weighted_seq =  graph_utils::BFS_forest(G, 1);
      updates = graph_utils::gen_random_weight_edges(non_weighted_seq, 1);
    }
    double time;
    std::string graph_name = file_name.substr(file_name.find_last_of('/')+1);
    graph_name = graph_name.substr(0, graph_name.length()-8);
    std::cout << "[ RUNNING " << graph_name << " MIS TREE BENCHMARK ]" << std::endl;
    output_csv << graph_name << "_mis" << ",";

    // RC Tree
    time = graph_benchmark::incremental_MSF_benchmark<TernarizedTree<RCTree<std::pair<int,Edge>>, std::pair<int,Edge>>>(updates, n);
    std::cout << "RCTree        : " << time << std::endl;
    output_csv << time << ",";
    // Topology Tree
    time = graph_benchmark::incremental_MSF_benchmark<TernarizedTree<TopologyTree<std::pair<int,Edge>, std::pair<int,Edge>>, std::pair<int,Edge>>>(updates, n);
    std::cout << "TopologyTree  : " << time << std::endl;
    output_csv << time << ",";
    // UFO Tree
    time = graph_benchmark::incremental_MSF_benchmark<UFOTree<std::pair<int,Edge>, std::pair<int,Edge>>>(updates, n);
    std::cout << "UFOTree       : " << time << std::endl;
    output_csv << time << ",";

    std::cout << std::endl;
    output_csv << "\n";
  }
}

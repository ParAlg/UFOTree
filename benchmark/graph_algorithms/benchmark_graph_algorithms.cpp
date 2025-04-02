#include "../dynamic_trees/benchmark.h"
#include "graph_benchmark.h"
#include "util.h"
#include "ternarized_tree.h"
#include "ufo_tree.h"
#include "topology_tree.h"
#include "rc_tree.h"
#include "parett/dynamic_trees/link_cut_tree/link_cut_tree.hpp"
#include "parett/dynamic_trees/euler_tour_tree/skip_list_ett.hpp"
#include "parett/dynamic_trees/euler_tour_tree/splay_tree_ett.hpp"
#include <fstream>


using namespace link_cut_tree;
using namespace skip_list_ett;

int main(int argc, char** argv) {

  srand(time(NULL));
  std::tuple<std::string, int> test_cases[] = {
    {"/ssd1/quinten/graphdata/com-youtube_sym.bin", 1},
    {"/ssd1/quinten/graphdata/as-skitter_sym.bin", 1},
    {"/ssd1/quinten/graphdata/enwiki_sym.bin", 1},
    {"/ssd1/quinten/graphdata/stackoverflow_sym.bin", 1},
    {"/ssd1/quinten/graphdata/RoadUSA_sym.bin", 1},
    {"/ssd1/quinten/graphdata/com-orkut_sym.bin", 1},
    {"/ssd1/quinten/graphdata/twitter_sym.bin", 1},
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
      auto unweighted_edges =  graph_utils::to_edges(G);
      updates = graph_utils::generate_random_weight_edges(unweighted_edges, 1);
    }
    double time;
    std::string graph_name = file_name.substr(file_name.find_last_of('/')+1);
    graph_name = graph_name.substr(0, graph_name.length()-8);
    std::cout << "[ RUNNING " << graph_name << " INCREMENTAL MSF BENCHMARK ]" << std::endl;
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
    // Link Cut Tree
    time = graph_benchmark::incremental_MSF_benchmark<LinkCutTreeInt>(updates, n);
    std::cout << "LinkCut Tree  : " << time << std::endl;
    output_csv << time << ",";

    std::cout << std::endl;
    output_csv << "\n";
  }
}

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

int main(int argc, char** argv) {
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

    // BFS trees
    for (auto test_case: test_cases) {
        std::string file_name = std::get<0>(test_case);
        int num_trials = std::get<1>(test_case);
        auto G = graph_utils::break_sym_graph_from_bin(file_name);
        auto n = G.size();
        std::vector<std::vector<Update>> update_sequences;
        for (int i = 0; i < num_trials; i++) {
            std::vector<Update> updates = graph_benchmark::stream_file_BFS_benchmark(G, rand());
            update_sequences.push_back(updates);
        }
        double time;
        std::string graph_name = file_name.substr(file_name.find_last_of('/')+1);
        graph_name = graph_name.substr(0, graph_name.length()-8);
        std::cout << "[ RUNNING " << graph_name << " BFS TREE UPDATE SPEED BENCHMARK ]" << std::endl;
        output_csv << graph_name << "_bfs" << ",";

        // RC Tree
        time = dynamic_tree_benchmark::get_update_speed<TernarizedTree<RCTree<int>, int>>(n, update_sequences);
        std::cout << "RCTree        : " << time << std::endl;
        output_csv << time << ",";
        // Topology Tree
        time = dynamic_tree_benchmark::get_update_speed<TernarizedTree<TopologyTree<int, empty_t>, empty_t>>(n, update_sequences);
        std::cout << "TopologyTree  : " << time << std::endl;
        output_csv << time << ",";
        // UFO Tree
        time = dynamic_tree_benchmark::get_update_speed<UFOTree<int, empty_t>>(n, update_sequences);
        std::cout << "UFOTree       : " << time << std::endl;
        output_csv << time << ",";
        // Euler Tour Tree
        time = dynamic_tree_benchmark::get_update_speed<EulerTourTree>(n, update_sequences);
        std::cout << "EulerTourTree : " << time << std::endl;
        output_csv << time << ",";

        std::cout << std::endl;
        output_csv << "\n";
    }

    // Random incremental spanning trees
    for (auto test_case: test_cases) {
        std::string file_name = std::get<0>(test_case);
        int num_trials = std::get<1>(test_case);
        auto G = graph_utils::break_sym_graph_from_bin(file_name);
        auto n = G.size();
        std::vector<std::vector<Update>> update_sequences;
        for (int i = 0; i < num_trials; i++) {
            std::vector<Update> updates = graph_benchmark::stream_file_RIS_benchmark(G, rand());
            update_sequences.push_back(updates);
        }
        double time;
        std::string graph_name = file_name.substr(file_name.find_last_of('/')+1);
        graph_name = graph_name.substr(0, graph_name.length()-8);
        std::cout << "[ RUNNING " << graph_name << " RIS TREE UPDATE SPEED BENCHMARK ]" << std::endl;
        output_csv << graph_name << "_ris" << ",";

        // RC Tree
        time = dynamic_tree_benchmark::get_update_speed<TernarizedTree<RCTree<int>, int>>(n, update_sequences);
        std::cout << "RCTree        : " << time << std::endl;
        output_csv << time << ",";
        // Topology Tree
        time = dynamic_tree_benchmark::get_update_speed<TernarizedTree<TopologyTree<int, empty_t>, empty_t>>(n, update_sequences);
        std::cout << "TopologyTree  : " << time << std::endl;
        output_csv << time << ",";
        // UFO Tree
        time = dynamic_tree_benchmark::get_update_speed<UFOTree<int, empty_t>>(n, update_sequences);
        std::cout << "UFOTree       : " << time << std::endl;
        output_csv << time << ",";
        // Euler Tour Tree
        time = dynamic_tree_benchmark::get_update_speed<EulerTourTree>(n, update_sequences);
        std::cout << "EulerTourTree : " << time << std::endl;
        output_csv << time << ",";

        std::cout << std::endl;
        output_csv << "\n";
    }

    // Incremental Minimum Spanning Tree
    
    Edge e; e.src = 0; e.dst = 1; std::pair<int,Edge> p(1, e);
    Edge e2; e.src = 1;e.dst = 2; std::pair<int,Edge> p2(3, e2);
    Edge e3; e.src = 0; e.dst = 2; std::pair<int,Edge> p3(2, e3);
    vector<std::pair<int, Edge>> v({p, p2, p3});
    graph_benchmark::incremental_MSF_benchmark<RCTree<std::pair<int, Edge>>>(v, 4);

    output_csv.close();
}

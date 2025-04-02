#include "../dynamic_trees/benchmark.h"
#include "graph_benchmark.h"
#include "util.h"
#include "ternarized_tree.h"
#include "ufo_tree.h"
#include "topology_tree.h"
#include "rc_tree.h"
#include "top_tree.h"
#include "parett/dynamic_trees/link_cut_tree/link_cut_tree.hpp"
#include "parett/dynamic_trees/euler_tour_tree/skip_list_ett.hpp"
#include "parett/dynamic_trees/euler_tour_tree/splay_tree_ett.hpp"
#include <fstream>


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

    // BFS trees
    output_csv << "Test Case,UFO Tree,Link Cut Tree,Splay Top Tree,Euler Tour Tree (Skip List),Euler Tour Tree (Splay Tree),Topology Tree,RC Tree";
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
        output_csv << "\n" << graph_name << "_bfs";

        // UFO Tree
        time = dynamic_tree_benchmark::get_update_speed<UFOTree<int, int>>(n, update_sequences);
        std::cout << "UFOTree       : " << time << std::endl;
        output_csv << "," << time;
        // Link Cut Tree
        time = dynamic_tree_benchmark::get_update_speed<link_cut_tree::LinkCutTreeInt>(n, update_sequences);
        std::cout << "LinkCutTree   : " << time << std::endl;
        output_csv << "," << time;
        // Top Tree
        time = dynamic_tree_benchmark::get_update_speed<TopTree<int>>(n, update_sequences);
        std::cout << "SplayTopTree  : " << time << std::endl;
        output_csv << "," << time;
        // Euler Tour Tree (Skip List)
        time = dynamic_tree_benchmark::get_update_speed<skip_list_ett::EulerTourTree>(n, update_sequences);
        std::cout << "SkipListETT   : " << time << std::endl;
        output_csv << "," << time;
        // Euler Tour Tree (Splay Tree)
        time = dynamic_tree_benchmark::get_update_speed<splay_tree_ett::EulerTourTree>(n, update_sequences);
        std::cout << "SplayTreeETT  : " << time << std::endl;
        output_csv << "," << time;
        // Topology Tree
        time = dynamic_tree_benchmark::get_update_speed<TernarizedTree<TopologyTree<int, int>, int>>(n, update_sequences);
        std::cout << "TopologyTree  : " << time << std::endl;
        output_csv << "," << time;
        // RC Tree
        time = dynamic_tree_benchmark::get_update_speed<TernarizedTree<RCTree<int>, int>>(n, update_sequences);
        std::cout << "RCTree        : " << time << std::endl;
        output_csv << "," << time;

        std::cout << std::endl;
    }

    // Random incremental spanning trees
    output_csv << "\n";
    output_csv << "Test Case,UFO Tree,Link Cut Tree,Splay Top Tree,Euler Tour Tree (Skip List),Euler Tour Tree (Splay Tree),Topology Tree,RC Tree";
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
        output_csv << "\n" << graph_name << "_ris";

        // UFO Tree
        time = dynamic_tree_benchmark::get_update_speed<UFOTree<int, int>>(n, update_sequences);
        std::cout << "UFOTree       : " << time << std::endl;
        output_csv << "," << time;
        // Link Cut Tree
        time = dynamic_tree_benchmark::get_update_speed<link_cut_tree::LinkCutTreeInt>(n, update_sequences);
        std::cout << "LinkCutTree   : " << time << std::endl;
        output_csv << "," << time;
        // Top Tree
        time = dynamic_tree_benchmark::get_update_speed<TopTree<int>>(n, update_sequences);
        std::cout << "SplayTopTree  : " << time << std::endl;
        output_csv << "," << time;
        // Euler Tour Tree (Skip List)
        time = dynamic_tree_benchmark::get_update_speed<skip_list_ett::EulerTourTree>(n, update_sequences);
        std::cout << "SkipListETT   : " << time << std::endl;
        output_csv << "," << time;
        // Euler Tour Tree (Splay Tree)
        time = dynamic_tree_benchmark::get_update_speed<splay_tree_ett::EulerTourTree>(n, update_sequences);
        std::cout << "SplayTreeETT  : " << time << std::endl;
        output_csv << "," << time;
        // Topology Tree
        time = dynamic_tree_benchmark::get_update_speed<TernarizedTree<TopologyTree<int, int>, int>>(n, update_sequences);
        std::cout << "TopologyTree  : " << time << std::endl;
        output_csv << "," << time;
        // RC Tree
        time = dynamic_tree_benchmark::get_update_speed<TernarizedTree<RCTree<int>, int>>(n, update_sequences);
        std::cout << "RCTree        : " << time << std::endl;
        output_csv << "," << time;

        std::cout << std::endl;
    }

    
    output_csv.close();
}

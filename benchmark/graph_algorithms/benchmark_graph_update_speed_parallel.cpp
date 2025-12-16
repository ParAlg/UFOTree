#include "../dynamic_trees/benchmark.h"
#include "graph_benchmark.h"
#include "parallel_benchmark.h"
#include "util.h"
#include "parallel_ufo_tree.h"
#include "parallel_topology_tree.h"
#include "parallel_topology_tree_ternarized.h"
#include "ParETT/euler_tour_tree.hpp"
#include "spaa_rc_tree.h"
#include "spaa_rc_tree_ternarized.h"
#include <fstream>


int main(int argc, char** argv) {
    // List of values of n to loop through and run all test cases
    vertex_t k = 10000;
    if (argc < 2) {
        std::cout << "Using default hard-coded value for k." << std::endl;
    } else if (argc == 2) {
        std::cout << "Using command line arguments for values of k." << std::endl;
        k = std::stoi(argv[1]);
    } else {
        std::cout << "Usage: ./parallel_graph_benchmark [k]" << std::endl;
    }
    long seed = 0;
    srand(seed);
    std::tuple<std::string, int> test_cases[] = {
        {"graphs/RoadUSA_sym.bin", 3},
        {"graphs/enwiki_sym.bin", 3},
        {"graphs/stackoverflow_sym.bin", 3},
        {"graphs/twitter_sym.bin", 3},
    };

    std::string filename = "../results/parallel_update_speed_graph.csv";
    std::ofstream output_csv;
    output_csv.open(filename);
    output_csv << "Test Case,Euler Tour Tree,UFO Tree,Topology Tree,Rake-Compress Tree,";

    // BFS trees
    for (auto test_case: test_cases) {
        std::string file_name = std::get<0>(test_case);
        int num_trials = std::get<1>(test_case);
        auto G = ufo::graph_utils::break_sym_graph_from_bin(file_name);
        auto n = G.size();
        std::vector<std::vector<UpdateBatch>> update_sequences(num_trials);
        std::vector<std::vector<UpdateBatchWithWeights>> weighted_update_sequences(num_trials);
        for (int i = 0; i < num_trials; i++) {
            std::vector<Update> updates = graph_benchmark::stream_file_BFS_benchmark(G, rand());
            update_sequences[i] = parallel_dynamic_tree_benchmark::convert_updates_to_batches(updates, k);
        }
        
        // Convert inputs to have dummy weights to fit Parallel RC-Tree Interface
        for (int i = 0; i < num_trials; i++) {
            weighted_update_sequences[i].resize(update_sequences[i].size());
            // Initialzing same update sequence with default val edge weights
            for (int j = 0; j < update_sequences[i].size(); j++) {
                auto updateBatch = update_sequences[i][j];
                auto edges = updateBatch.edges;
                weighted_update_sequences[i][j].type = updateBatch.type;
                if (updateBatch.type == INSERT) {
                    for (int k = 0; k < edges.size(); k++) {
                        weighted_update_sequences[i][j].insert_edges.push_back(std::make_tuple(edges[k].first, edges[k].second, 0));
                    }
                } else {
                    weighted_update_sequences[i][j].delete_edges = edges;
                }
            }
        }

        double time;
        std::string graph_name = file_name.substr(file_name.find_last_of('/')+1);
        graph_name = graph_name.substr(0, graph_name.length()-8);
        std::cout << "[ RUNNING " << graph_name << " BFS TREE UPDATE SPEED BENCHMARK ]" << std::endl;
        output_csv << "\n" << graph_name << "_bfs,";

        // Euler Tour Tree
        time = parallel_dynamic_tree_benchmark::get_update_speed<parallel_euler_tour_tree::EulerTourTree<int>>(n, k, update_sequences);
        std::cout << "EulerTourTree : " << time << std::endl;
        output_csv << time << ",";
        // UFO Tree
        time = parallel_dynamic_tree_benchmark::get_update_speed<ParallelUFOTree<>>(n, k, update_sequences);
        std::cout << "UFOTree       : " << time << std::endl;
        output_csv << time << ",";
        // Topology Tree
        time = parallel_dynamic_tree_benchmark::get_update_speed_with_rand_edge_weights<ParallelTopologyTreeTernarized<int>>(n, k, weighted_update_sequences);
        std::cout << "TopologyTree  : " << time << std::endl;
        output_csv << time << ",";
        // RC Tree
        time = parallel_dynamic_tree_benchmark::get_update_speed_with_rand_edge_weights<ParallelRCTreeTernarized<int>>(n, k, weighted_update_sequences);
        std::cout << "RCTree        : " << time << std::endl;
        output_csv << time << ",";
        std::cout << std::endl;
    }

    // Random incremental spanning trees
    for (auto test_case: test_cases) {
        std::string file_name = std::get<0>(test_case);
        int num_trials = std::get<1>(test_case);
        auto G = ufo::graph_utils::break_sym_graph_from_bin(file_name);
        auto n = G.size();
        std::vector<std::vector<UpdateBatch>> update_sequences(num_trials);
        std::vector<std::vector<UpdateBatchWithWeights>> weighted_update_sequences(num_trials);
        for (int i = 0; i < num_trials; i++) {
            std::vector<Update> updates = graph_benchmark::stream_file_RIS_benchmark(G, rand());
            update_sequences[i] = parallel_dynamic_tree_benchmark::convert_updates_to_batches(updates, k);
        }
        
        // Convert inputs to have dummy weights to fit Parallel RC-Tree Interface
        for (int i = 0; i < num_trials; i++) {
            weighted_update_sequences[i].resize(update_sequences[i].size());
            // Initialzing same update sequence with default val edge weights
            for (int j = 0; j < update_sequences[i].size(); j++) {
                auto updateBatch = update_sequences[i][j];
                auto edges = updateBatch.edges;
                weighted_update_sequences[i][j].type = updateBatch.type;
                if (updateBatch.type == INSERT) {
                    for (int k = 0; k < edges.size(); k++) {
                        weighted_update_sequences[i][j].insert_edges.push_back(std::make_tuple(edges[k].first, edges[k].second, 0));
                    }
                } else {
                    weighted_update_sequences[i][j].delete_edges = edges;
                }
            }
        }

        double time;
        std::string graph_name = file_name.substr(file_name.find_last_of('/')+1);
        graph_name = graph_name.substr(0, graph_name.length()-8);
        std::cout << "[ RUNNING " << graph_name << " RIS TREE UPDATE SPEED BENCHMARK ]" << std::endl;
        output_csv << "\n" << graph_name << "_ris,";

        // Euler Tour Tree
        time = parallel_dynamic_tree_benchmark::get_update_speed<parallel_euler_tour_tree::EulerTourTree<int>>(n, k, update_sequences);
        std::cout << "EulerTourTree : " << time << std::endl;
        output_csv << time << ",";
        // UFO Tree
        time = parallel_dynamic_tree_benchmark::get_update_speed<ParallelUFOTree<>>(n, k, update_sequences);
        std::cout << "UFOTree       : " << time << std::endl;
        output_csv << time << ",";
        // Topology Tree
        time = parallel_dynamic_tree_benchmark::get_update_speed_with_rand_edge_weights<ParallelTopologyTreeTernarized<int>>(n, k, weighted_update_sequences);
        std::cout << "TopologyTree  : " << time << std::endl;
        output_csv << time << ",";
        // RC Tree
        time = parallel_dynamic_tree_benchmark::get_update_speed_with_rand_edge_weights<ParallelRCTreeTernarized<int>>(n, k, weighted_update_sequences);
        std::cout << "RCTree        : " << time << std::endl;
        output_csv << time << ",";
        std::cout << std::endl;

        std::cout << std::endl;
    }

    
    output_csv.close();
}

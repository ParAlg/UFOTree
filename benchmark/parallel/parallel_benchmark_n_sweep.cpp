#include "../dynamic_trees/benchmark.h"
#include "parallel_benchmark.h"
#include "util.h"
#include "parallel_ufo_tree.h"
#include "parallel_topology_tree.h"
#include "parallel_topology_tree_ternarized.h"
#include "ParETT/euler_tour_tree.hpp"
#include "spaa_rc_tree.h"
#include "spaa_rc_tree_ternarized.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <functional>

using namespace ufo;

int main(int argc, char** argv) {
    // List of values of n to loop
    std::vector<vertex_t> n_list = {1000001, 5000000, 10000000, 50000000, 100000000, 500000000, 1000000000};
    // std::vector<vertex_t> n_list = {10000, 50000, 100000};
    vertex_t k = 1E6;
    long seed = 0;
    srand(seed);

    /* Define the specific test cases requested for the columns:
       1. Path
       2. Binary
       3. 64ary
       4. Star
    */
    std::tuple<std::string, std::function<std::vector<Update>(vertex_t, long)>, bool, int> test_cases[] = {
        {"Path", dynamic_tree_benchmark::linked_list_benchmark, false, 1},
        {"Binary", dynamic_tree_benchmark::binary_tree_benchmark, false, 1},
        {"64ary", dynamic_tree_benchmark::k_ary_tree_benchmark, true, 1},
        {"Star", dynamic_tree_benchmark::star_benchmark, true, 1},
    };

    // Open a single shared CSV file
    std::string filename = "../results/n_sweep_parallel_update.csv";
    std::ofstream output_csv;
    output_csv.open(filename);

    if (!output_csv.is_open()) {
        std::cerr << "Error: Could not open output file " << filename << std::endl;
        return 1;
    }

    // Write the header row
    output_csv << "n,Path,Binary,64ary,Star\n";

    // Outer loop: Iterate over input sizes (n)
    for (vertex_t n : n_list) {
        std::cout << "[ Processing n=" << n << " ]" << std::endl;
        
        // Start the CSV row with the current n
        output_csv << n;

        // Inner loop: Run each specific benchmark type for this n
        for (auto& test_case : test_cases) {
            std::string test_case_name = std::get<0>(test_case);
            auto update_generator = std::get<1>(test_case);
            bool ternarize = std::get<2>(test_case);
            int num_trials = std::get<3>(test_case);

            std::vector<std::vector<UpdateBatch>> update_sequences(num_trials);
            std::vector<std::vector<UpdateBatchWithWeights>> weighted_update_sequences(num_trials);

            // Generate updates for this specific n and structure
            for (int i = 0; i < num_trials; i++) {
                auto update_sequence = update_generator(n, rand());
                update_sequences[i] = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);
                
                // Keep weighted logic in case RC trees are re-enabled later
                weighted_update_sequences[i].resize(update_sequences[i].size());
                for(size_t j = 0; j < update_sequences[i].size(); j++){
                    auto updateBatch = update_sequences[i][j];
                    auto edges = updateBatch.edges;
                    weighted_update_sequences[i][j].type = updateBatch.type;
                    if (updateBatch.type == INSERT){
                        for (size_t k_edge = 0; k_edge < edges.size(); k_edge++){
                            weighted_update_sequences[i][j].insert_edges.push_back(std::make_tuple(edges[k_edge].first, edges[k_edge].second, 0));
                        }
                    } else {
                        weighted_update_sequences[i][j].delete_edges = edges;
                    }
                }
            }

            // Run UFO Tree Benchmark
            std::cout << "  Running " << test_case_name << "... ";
            double time = parallel_dynamic_tree_benchmark::get_update_speed<ParallelUFOTree<>>(n, k, update_sequences);
            std::cout << time << "s" << std::endl;

            // Append time to the current CSV row
            output_csv << "," << time;
        }

        // End the CSV row
        output_csv << "\n";
    }

    output_csv.close();
    std::cout << "Results saved to " << filename << std::endl;

    return 0;
}
#include "../dynamic_trees/benchmark.h"
#include "parallel_benchmark.h"
#include "util.h"
#include "parallel_ufo_tree.h"
#include "parallel_topology_tree.h"
#include "ParETT/euler_tour_tree.hpp"
#include "spaa_rc_tree.h"
#include <fstream>


using namespace dgbs;


parlay::internal::timer timer0("");
parlay::internal::timer timer1("");
parlay::internal::timer timer2("");
parlay::internal::timer timer3("");
parlay::internal::timer timer4("");
parlay::internal::timer timer5("");

parlay::internal::timer subtimer1("");
parlay::internal::timer subtimer2("");
parlay::internal::timer subtimer3("");
parlay::internal::timer subtimer4("");

int main(int argc, char** argv) {
    // List of values of n to loop through and run all test cases
    vertex_t n = 100000;
    vertex_t k = 10000;
    if (argc < 2) {
        std::cout << "Using default hard-coded list for values of n and k." << std::endl;
    } else if (argc == 3) {
        std::cout << "Using command line arguments for values of n and k." << std::endl;
        n = std::stoi(argv[1]);
        k = std::stoi(argv[2]);
    } else {
        std::cout << "Usage: ./parallel_benchmark [n] [k]" << std::endl;
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

    std::string filename = "../results/parellel_update_speed_" + std::to_string(n) + "_" + std::to_string(k) + ".csv";
    std::ofstream output_csv;
    output_csv.open(filename);
    output_csv << "Test Case,Euler Tour Tree,UFO Tree,Topology Tree,\n";

    for (auto test_case : test_cases) {
        std::string test_case_name = std::get<0>(test_case);
        auto update_generator = std::get<1>(test_case);
        bool ternarize = std::get<2>(test_case);
        int num_trials = std::get<3>(test_case);
        std::vector<std::vector<UpdateBatch>> update_sequences(num_trials);
        std::vector<std::vector<UpdateBatchWithWeights>> weighted_update_sequences(num_trials);

        //Generate weighted Update Sequence for RC trees (only in the case of insertion)
        for (int i = 0; i < num_trials; i++) {
            auto update_sequence = update_generator(n, rand());
            update_sequences[i] = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);
            weighted_update_sequences[i].resize(update_sequences[i].size());
            // Initialzing same update sequence with default val edge weights
            for(int j = 0; j < update_sequences[i].size(); j++){
                auto updateBatch = update_sequences[i][j];
                auto edges = updateBatch.edges;
                weighted_update_sequences[i][j].type = updateBatch.type;
                if (updateBatch.type == INSERT){
                for (int k = 0; k < edges.size(); k++){
                    weighted_update_sequences[i][j].insert_edges.push_back(make_tuple(edges[k].first, edges[k].second, 0));
                }
                } else {
                    weighted_update_sequences[i][j].delete_edges = edges;
                }
            }
        }
        double time;
        std::cout << "[ RUNNING " << test_case_name << " PARALLEL UPDATE SPEED BENCHMARK WITH n=" << n << ", k=" << k << " ]" << std::endl;
        output_csv << test_case_name << ",";

        // UFO Tree
        time = parallel_dynamic_tree_benchmark::get_update_speed<ParallelUFOTree<>>(n, k, update_sequences);
        std::cout << "UFOTree       : " << time << std::endl;
        output_csv << time << ",";
        // Euler Tour Tree
        // time = parallel_dynamic_tree_benchmark::get_update_speed<parallel_euler_tour_tree::EulerTourTree<int>>(n, k, update_sequences);
        // std::cout << "EulerTourTree : " << time << std::endl;
        // output_csv << time << ",";
        // Topology Tree
        // if (!ternarize) {
        //     time = parallel_dynamic_tree_benchmark::get_update_speed<ParallelTopologyTree<int>>(n, k, update_sequences);
        //     std::cout << "TopologyTree  : " << time << std::endl;
        //     output_csv << time << ",";
        // }
        // RC Tree
        // if (!ternarize) {
        //     time = parallel_dynamic_tree_benchmark::get_update_speed_with_rand_edge_weights<ParallelRCTree<int>>(n,k,weighted_update_sequences);
        //     std::cout << "RCTree        : " << time << std::endl;
        //     output_csv << time << ",";
        // }
        std::cout << std::endl;
        output_csv << "\n";
    }

    output_csv.close();

    std::cout << "TIMER 0: " << timer0.total_time() << std::endl;
    std::cout << "TIMER 1: " << timer1.total_time() << std::endl;
    std::cout << "TIMER 2: " << timer2.total_time() << std::endl;
    std::cout << "TIMER 3: " << timer3.total_time() << std::endl;
    std::cout << "TIMER 4: " << timer4.total_time() << std::endl;
    std::cout << "TIMER 5: " << timer5.total_time() << std::endl;

    std::cout << "SUBTIMER 1: " << subtimer1.total_time() << std::endl;
    std::cout << "SUBTIMER 2: " << subtimer2.total_time() << std::endl;
    std::cout << "SUBTIMER 3: " << subtimer3.total_time() << std::endl;
    std::cout << "SUBTIMER 4: " << subtimer4.total_time() << std::endl;
}

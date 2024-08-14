#include "benchmark.h"
#include "../include/ufo_tree_abseil.h"
#include "../include/topology_tree.h"
#include "../include/rc_tree.h"
#include "../baselines/dynamic_trees/euler_tour_tree/include/skip_list_ett.hpp"


int main(int argc, char** argv) {
    // Run each test case using each data structure
    std::string test_case;
    // vertex_t n_list[] = {1000000, 10000000};

    for (vertex_t n : n_list) {

        test_case = "linked-list";
        std::cout << "[ RUNNING " << test_case << " BENCHMARK WITH n=" << n << " ]" << std::endl;
        std::cout << "RCTree        ";
        dynamic_tree_benchmark::linked_list_benchmark<RCTree<int>>(n);
        std::cout << "TopologyTree  ";
        dynamic_tree_benchmark::linked_list_benchmark<TopologyTree<int>>(n);
        std::cout << "UFOTree       ";
        dynamic_tree_benchmark::linked_list_benchmark<UFOTree<int>>(n);
        std::cout << "EulerTourTree ";
        dynamic_tree_benchmark::linked_list_benchmark<skip_list_ett::EulerTourTree>(n);
        std::cout << std::endl;

        test_case = "binary-tree";
        std::cout << "[ RUNNING " << test_case << " BENCHMARK WITH n=" << n << " ]" << std::endl;
        std::cout << "RCTree        ";
        dynamic_tree_benchmark::binary_tree_benchmark<RCTree<int>>(n);
        std::cout << "TopologyTree  ";
        dynamic_tree_benchmark::binary_tree_benchmark<TopologyTree<int>>(n);
        std::cout << "UFOTree       ";
        dynamic_tree_benchmark::binary_tree_benchmark<UFOTree<int>>(n);
        std::cout << "EulerTourTree ";
        dynamic_tree_benchmark::binary_tree_benchmark<skip_list_ett::EulerTourTree>(n);
        std::cout << std::endl;

        test_case = "64ary-tree";
        std::cout << "[ RUNNING " << test_case << " BENCHMARK WITH n=" << n << " ]" << std::endl;
        std::cout << "UFOTree       ";
        dynamic_tree_benchmark::k_ary_tree_benchmark<UFOTree<int>>(n);
        std::cout << "EulerTourTree ";
        dynamic_tree_benchmark::k_ary_tree_benchmark<skip_list_ett::EulerTourTree>(n);
        std::cout << std::endl;

        test_case = "star";
        std::cout << "[ RUNNING " << test_case << " BENCHMARK WITH n=" << n << " ]" << std::endl;
        std::cout << "UFOTree       ";
        dynamic_tree_benchmark::star_benchmark<UFOTree<int>>(n);
        std::cout << "EulerTourTree ";
        dynamic_tree_benchmark::star_benchmark<skip_list_ett::EulerTourTree>(n);
        std::cout << std::endl;

        test_case = "random-degree-3";
        std::cout << "[ RUNNING " << test_case << " BENCHMARK WITH n=" << n << " ]" << std::endl;
        std::cout << "RCTree        ";
        dynamic_tree_benchmark::random_degree3_benchmark<RCTree<int>>(n);
        std::cout << "TopologyTree  ";
        dynamic_tree_benchmark::random_degree3_benchmark<TopologyTree<int>>(n);
        std::cout << "UFOTree       ";
        dynamic_tree_benchmark::random_degree3_benchmark<UFOTree<int>>(n);
        std::cout << "EulerTourTree ";
        dynamic_tree_benchmark::random_degree3_benchmark<skip_list_ett::EulerTourTree>(n);
        std::cout << std::endl;

        test_case = "random-unbounded";
        std::cout << "[ RUNNING " << test_case << " BENCHMARK WITH n=" << n << " ]" << std::endl;
        std::cout << "UFOTree       ";
        dynamic_tree_benchmark::random_unbounded_benchmark<UFOTree<int>>(n);
        std::cout << "EulerTourTree ";
        dynamic_tree_benchmark::random_unbounded_benchmark<skip_list_ett::EulerTourTree>(n);
        std::cout << std::endl;

        /*test_case = "preferential-attachment";
        std::cout << "[ RUNNING " << test_case << " BENCHMARK WITH n=" << n << " ]" << std::endl;
        std::cout << "UFOTree       ";
        dynamic_tree_benchmark::preferential_attachment_benchmark<UFOTree<int>>(n);
        std::cout << "EulerTourTree ";
        dynamic_tree_benchmark::preferential_attachment_benchmark<skip_list_ett::EulerTourTree>(n);
        std::cout << std::endl;*/
    }
}

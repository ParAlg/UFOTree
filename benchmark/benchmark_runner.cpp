#include "benchmark.h"
#include "../include/ufo_tree.h"
#include "../include/topology_tree.h"
#include "../include/rc_tree.h"
#include "../baselines/dynamic_trees/euler_tour_tree/include/skip_list_ett.hpp"


int main(int argc, char** argv) {
    // Run each test case using each data structure
    auto test_case = "";
    auto data_structure = "";
    vertex_t n = 1000000;


    test_case = "incremental-linked-list";

    data_structure = "UFOTree";
    std::cout << std::endl << "RUNNING: [ " << data_structure << " " << test_case << " n=" << n << " ]" << std::endl;
    dynamic_tree_benchmark::incremental_linked_list_benchmark<UFOTree<int>>(n);

    data_structure = "TopologyTree";
    std::cout << std::endl << "RUNNING: [ " << data_structure << " " << test_case << " n=" << n << " ]" << std::endl;
    dynamic_tree_benchmark::incremental_linked_list_benchmark<TopologyTree<int>>(n);

    data_structure = "EulerTourTree";
    std::cout << std::endl << "RUNNING: [ " << data_structure << " " << test_case << " n=" << n << " ]" << std::endl;
    dynamic_tree_benchmark::incremental_linked_list_benchmark<skip_list_ett::EulerTourTree>(n);


    test_case = "random-degree3";

    data_structure = "UFOTree";
    std::cout << std::endl << "RUNNING: [ " << data_structure << " " << test_case << " n=" << n << " ]" << std::endl;
    dynamic_tree_benchmark::random_degree3_benchmark<UFOTree<int>>(n);

    data_structure = "TopologyTree";
    std::cout << std::endl << "RUNNING: [ " << data_structure << " " << test_case << " n=" << n << " ]" << std::endl;
    dynamic_tree_benchmark::random_degree3_benchmark<TopologyTree<int>>(n);

    data_structure = "EulerTourTree";
    std::cout << std::endl << "RUNNING: [ " << data_structure << " " << test_case << " n=" << n << " ]" << std::endl;
    dynamic_tree_benchmark::random_degree3_benchmark<skip_list_ett::EulerTourTree>(n);


    test_case = "random-unbounded";

    data_structure = "UFOTree";
    std::cout << std::endl << "RUNNING: [ " << data_structure << " " << test_case << " n=" << n << " ]" << std::endl;
    dynamic_tree_benchmark::random_unbounded_benchmark<UFOTree<int>>(n);

    data_structure = "EulerTourTree";
    std::cout << std::endl << "RUNNING: [ " << data_structure << " " << test_case << " n=" << n << " ]" << std::endl;
    dynamic_tree_benchmark::random_unbounded_benchmark<skip_list_ett::EulerTourTree>(n);
}

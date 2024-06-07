#include "benchmark.h"
#include "../include/ufo_tree.h"
#include "../include/topology_tree.h"
#include "../include/rc_tree.h"
#include "../baselines/dynamic_trees/euler_tour_tree/include/skip_list_ett.hpp"


int main(int argc, char** argv) {
    vertex_t n = 1000000;
    std::cout << std::endl << "RUNNING: [ " << "UFOTree" << "\t " << "incremental-linked-list" << "\t n=" << n << " ]" << std::endl;
    dynamic_tree_benchmark::run_sequential_benchmark<UFOTree<int>>(n);
    std::cout << std::endl << "RUNNING: [ " << "TopologyTree" << "\t " << "incremental-linked-list" << "\t n=" << n << " ]" << std::endl;
    dynamic_tree_benchmark::run_sequential_benchmark<TopologyTree<int>>(n);
}

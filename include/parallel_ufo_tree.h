#pragma once
#include "types.h"
#include "util.h"
#include "bridge.h"
#include "hash_bag.h"
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include "sequential_forest.h"


// #define PRINT_DEBUG_INFO

using namespace parlay;

namespace dgbs {

struct ParallelUFOCluster {
    // Topology cluster data
    std::unordered_set<ParallelUFOCluster*> neighbors;
    ParallelUFOCluster* parent;

    // Constructor
    ParallelUFOCluster(){};
    // Helper functions
    int get_degree();
    bool contracts();
    bool contains_neighbor(ParallelUFOCluster* c);
    void insert_neighbor(ParallelUFOCluster* c);
    void remove_neighbor(ParallelUFOCluster* c);
    ParallelUFOCluster* get_root();
};

template <typename aug_t = empty_t>
class ParallelUFOTree {
public:
    // UFO tree interface
    ParallelUFOTree(vertex_t n, vertex_t k);
    void batch_link(parlay::sequence<std::pair<int, int>>& links);
    void batch_cut(parlay::sequence<std::pair<int, int>>& cuts);
    bool connected(vertex_t u, vertex_t v);
    // Testing helpers
    bool is_valid();
    int get_height(vertex_t v);
    void print_tree();

private:
    // Class data and parameters
    vertex_t n;
    ParallelUFOCluster* leaves;
    // Main batch update functions
    void recluster_tree();
};

template <typename aug_t>
ParallelUFOTree<aug_t>::ParallelUFOTree(vertex_t n, vertex_t k) : n(n) {
    leaves = new ParallelUFOCluster[n];
    for (uint32_t i = 0; i < n; ++i)
        leaves[i] = ParallelUFOCluster();
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_link(parlay::sequence<std::pair<int, int>>& links) {
    
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_cut(parlay::sequence<std::pair<int, int>>& cuts) {
    
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_tree() {
    
}

template <typename aug_t>
bool ParallelUFOTree<aug_t>::connected(vertex_t u, vertex_t v) {

}

}

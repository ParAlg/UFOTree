#pragma once

#include "types.h"
#include "util.h"
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include "sequential_forest.h"


using namespace parlay;

template <typename aug_t>
class ParallelUFOTree {
public:
    using ForestStructureType = SequentialForest;
    // UFO tree interface
    ParallelUFOTree(vertex_t n, vertex_t k, QueryType q = PATH,
        std::function<aug_t(aug_t, aug_t)> f = [](aug_t x, aug_t y) -> aug_t {
        return x + y; }, aug_t id = 0, aug_t dval = 0);
    void batch_link(Edge* links, int len);
    void batch_cut(Edge* cuts, int len);
    bool connected(vertex_t u, vertex_t v);
    // Testing helpers
    bool is_valid();
    int get_height(vertex_t v);
    void print_tree();

private:
    // Class data and parameters
    sequence<ForestStructureType> levels;
    vertex_t n;
    QueryType query_type;
    std::function<aug_t(aug_t, aug_t)> f;
    aug_t identity;
    aug_t default_value;
    // Helper functions
    void recluster_tree();
};

template <typename aug_t>
ParallelUFOTree<aug_t>::ParallelUFOTree(vertex_t n, vertex_t k, QueryType q,
std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d)
    : n(n), query_type(q), f(f), identity(id), default_value(d) {
    ForestStructureType F_0;
    sequence<vertex_t> leaves = tabulate(n, [&] (size_t i) { return (vertex_t) i; });
    F_0.insert_vertices(leaves);
    levels.push_back(F_0);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_link(Edge* links, int len) {

}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_cut(Edge* cuts, int len) {
    
}

template <typename aug_t>
bool ParallelUFOTree<aug_t>::connected(vertex_t u, vertex_t v) {
    vertex_t root_u = u;
    int level = 0;
    while (levels[level].get_parent(root_u) != -1)
        root_u = levels[level++].get_parent(root_u);
    vertex_t root_v = v;
    level = 0;
    while (levels[level].get_parent(root_v) != -1)
        root_v = levels[level++].get_parent(root_v);
    return root_u == root_v;
}

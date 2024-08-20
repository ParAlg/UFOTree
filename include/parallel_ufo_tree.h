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
    void batch_link(sequence<Edge>& links);
    void batch_cut(sequence<Edge>& cuts);
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
    void update_tree(sequence<Edge>& updates, bool deletion);
    void recluster_level(int level, bool deletion, sequence<vertex_t>& R, sequence<vertex_t>& D, sequence<Edge>& U);
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
void ParallelUFOTree<aug_t>::batch_link(sequence<Edge>& links) {
    update_tree(links, false);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_cut(sequence<Edge>& cuts) {
    update_tree(cuts, true);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::update_tree(sequence<Edge>& updates, bool deletion) {
    // Determine the intial set of root clusters and deletion clusters
    sequence<vertex_t> R = levels[0].get_endpoints(updates);
    sequence<vertex_t> D = levels[0].get_parents(R);
    // Remove parents of low degree level 0 root clusters
    sequence<vertex_t> low_deg = filter(R, [&] (auto v) { return levels[0].get_degree(v) < 2; });
    levels[0].unset_parents(low_deg);
    sequence<vertex_t> low_deg_parents = parlay::map(low_deg, [&] (auto v) { return levels[0].get_parent(v); });
    if (levels.size() > 1) levels[1].subtract_children(low_deg_parents);
    // Perform the updates at level 0 and determine the next level updates
    sequence<Edge> U;
    if (!deletion) {
        levels[0].insert_edges(updates);
        U = updates;
    } else {
        levels[0].delete_edges(updates);
        U = levels[0].filter_edges(updates);
    }
    // Start the level-by-level update procedure
    recluster_level(0, deletion, R, D, U);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_level(int level, bool deletion, sequence<vertex_t>& R, sequence<vertex_t>& D, sequence<Edge>& U) {

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

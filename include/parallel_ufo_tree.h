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
    levels.emplace_back();
    sequence<vertex_t> leaves = tabulate(n, [&] (size_t i) { return (vertex_t) i; });
    levels[0].insert_vertices(leaves);
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
    sequence<vertex_t> low_deg_parents = parlay::map(low_deg, [&] (auto v) { return levels[0].get_parent(v); });
    levels[0].unset_parents(low_deg);
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
    levels.pop_back();
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_level(int level, bool deletion, sequence<vertex_t>& R, sequence<vertex_t>& D, sequence<Edge>& U) {
    // Delete the edges from our initial updates that are still in level i+1
    if (levels.size() < level) levels.emplace_back();
    if (deletion) levels[level+1].delete_edges(U);

    // Determine the set of clusters that will be deleted in level i+1
    sequence<vertex_t> del = filter(D, [&] (auto v) {
        bool low_degree = levels[level+1].get_degree(v) < 3;
        bool low_fanout = levels[level+1].get_child_count(v) < 3;
        return (low_degree && low_fanout);
    });
    // Collect the new level i root clusters that will be formed by deleting some level i+1 clusters and the parents
    // levels[level].getSiblings(R);
    sequence<vertex_t> next_D = levels[level+1].get_parents(D);
    // Do the deletion of clusters
    levels[level+1].delete_vertices(del);

    // Remove the parents of low degree root clusters
    sequence<vertex_t> low_deg = filter(R, [&] (auto v) { return levels[level].get_degree(v) < 2; });
    sequence<vertex_t> low_deg_parents = parlay::map(low_deg, [&] (auto v) { return levels[level].get_parent(v); });
    levels[level].unset_parents(low_deg);
    levels[level+1].subtract_children(low_deg_parents);

    // Determine the new combinations of the root clusters at this level
    // levels[level].MaximalCombination(R);

    // Prepare the inputs to the next level
    sequence<vertex_t> next_R = levels[level].get_parents(R);
    sequence<vertex_t> additional_D = levels[level+1].get_parents(next_R);
    next_D = remove_duplicates(append(next_D, additional_D));
    sequence<Edge> next_U = levels[level+1].filter_edges(U);
    // Remove the parents of low degree clusters in the next root clusters
    low_deg = filter(next_R, [&] (auto v) { return levels[level+1].get_degree(v) < 2; });
    low_deg_parents = parlay::map(low_deg, [&] (auto v) { return levels[level+1].get_parent(v); });
    levels[level+1].unset_parents(low_deg);
    if (levels.size() > level+2) levels[level+2].subtract_children(low_deg_parents);
    // Check if the update is complete or recursively recluster the next level
    if (next_R.empty() && next_D.empty()) return;
    recluster_level(level+1, deletion, next_R, next_D, next_U);
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

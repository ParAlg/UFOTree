#pragma once

#include "types.h"
#include "util.h"
#include "bridge.h"
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
    sequence<ForestStructureType> forests;
    vertex_t n;
    QueryType query_type;
    std::function<aug_t(aug_t, aug_t)> f;
    aug_t identity;
    aug_t default_value;
    // Helper functions
    void update_tree(sequence<Edge>& updates, bool deletion);
    void recluster_level(int level, bool deletion, sequence<vertex_t>& R, sequence<vertex_t>& D, sequence<Edge>& U);
    sequence<std::pair<vertex_t,vertex_t>> maximal_combination(int level, sequence<vertex_t>& R);
};

template <typename aug_t>
ParallelUFOTree<aug_t>::ParallelUFOTree(vertex_t n, vertex_t k, QueryType q,
std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d)
    : n(n), query_type(q), f(f), identity(id), default_value(d) {
    forests.emplace_back();
    sequence<vertex_t> leaves = tabulate(n, [&] (size_t i) { return (vertex_t) i; });
    forests[0].insert_vertices(leaves);
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
    sequence<vertex_t> R = forests[0].get_endpoints(updates);
    sequence<vertex_t> D = forests[0].get_parents(R);
    // Remove parents of low degree level 0 root clusters
    sequence<vertex_t> low_deg = filter(R, [&] (auto v) { return forests[0].get_degree(v) < 2; });
    sequence<vertex_t> low_deg_parents = parlay::map(low_deg, [&] (auto v) { return forests[0].get_parent(v); });
    forests[0].unset_parents(low_deg);
    if (forests.size() > 1) forests[1].subtract_children(low_deg_parents);
    // Perform the updates at level 0 and determine the next level updates
    sequence<Edge> U;
    if (!deletion) {
        forests[0].insert_edges(updates);
        U = updates;
    } else {
        forests[0].delete_edges(updates);
        U = forests[0].filter_edges(updates);
    }
    // Start the level-by-level update procedure
    recluster_level(0, deletion, R, D, U);
    forests.pop_back();
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_level(int level, bool deletion, sequence<vertex_t>& R, sequence<vertex_t>& D, sequence<Edge>& U) {
    // Delete the edges from our initial updates that are still in level i+1
    if (forests.size() < level) forests.emplace_back();
    if (deletion) forests[level+1].delete_edges(U);

    // Determine the set of clusters that will be deleted in level i+1
    sequence<vertex_t> del = filter(D, [&] (auto v) {
        bool low_degree = forests[level+1].get_degree(v) < 3;
        bool low_fanout = forests[level+1].get_child_count(v) < 3;
        return (low_degree && low_fanout);
    });
    // Collect the new level i root clusters that will be formed by deleting some level i+1 clusters and the parents
    // forests[level].getSiblings(R);
    sequence<vertex_t> next_D = forests[level+1].get_parents(D);
    // Do the deletion of clusters
    forests[level+1].delete_vertices(del);

    // Remove the parents of low degree root clusters
    sequence<vertex_t> low_deg = filter(R, [&] (auto v) { return forests[level].get_degree(v) < 2; });
    sequence<vertex_t> low_deg_parents = parlay::map(low_deg, [&] (auto v) { return forests[level].get_parent(v); });
    forests[level].unset_parents(low_deg);
    forests[level+1].subtract_children(low_deg_parents);

    // Determine the new combinations of the root clusters at this level
    auto combinations = maximal_combination(level, R);

    // Prepare the inputs to the next level
    sequence<vertex_t> next_R = forests[level].get_parents(R);
    sequence<vertex_t> additional_D = forests[level+1].get_parents(next_R);
    next_D = remove_duplicates(append(next_D, additional_D));
    sequence<Edge> next_U = forests[level+1].filter_edges(U);
    // Remove the parents of low degree clusters in the next root clusters
    low_deg = filter(next_R, [&] (auto v) { return forests[level+1].get_degree(v) < 2; });
    low_deg_parents = parlay::map(low_deg, [&] (auto v) { return forests[level+1].get_parent(v); });
    forests[level+1].unset_parents(low_deg);
    if (forests.size() > level+2) forests[level+2].subtract_children(low_deg_parents);
    // Check if the update is complete or recursively recluster the next level
    if (next_R.empty() && next_D.empty()) return;
    recluster_level(level+1, deletion, next_R, next_D, next_U);
}

template <typename aug_t>
sequence<std::pair<vertex_t,vertex_t>> ParallelUFOTree<aug_t>::maximal_combination(int level, sequence<vertex_t>& R) {
    // Perform clustering of the root clusters
    parallel_for(0, R.size(), [&](size_t i) {
        auto cluster = R[i];
        if (forests[level].get_parent(cluster)) return;
        // Combine deg 3+ root clusters with deg 1 root or non-root clusters
        if (forests[level].get_degree(cluster) >= 3) {
            // If it became high degree, any non-contracting deg 1 clusters must be in the first 2 neighbors
            for (auto neighbor : cluster->neighbors) {
                if (forests[level].get_degree(neighbor) == 1) {
                    forests[level].try_set_partner(cluster, neighbor);
                    forests[level].try_set_partner(neighbor, cluster);
                    forests[level].unset_parent(neighbor);
                }
            }
        } else if (forests[level].get_degree(cluster) == 2) {
            // Only local maxima in priority with respect to deg 2 clusters will act
            bool local_max = true;
            for (auto neighbor : cluster->neighbors) {
                if (!neighbor->parent && forests[level].get_degree(neighbor) == 2
                && neighbor->priority >= cluster->priority) {
                    local_max = false;
                }
            }
            if (!local_max) return;
            // Travel left/right and pair clusters until a deg 3, deg 1, non-root,
            // or combined cluster found
            for (auto direction : {0, 1}) {
                auto curr = cluster;
                auto next = cluster->get_neighbor(direction);
                if (curr->partner) {
                    curr = cluster->get_neighbor(direction);
                    next = curr->get_other_neighbor(cluster);
                }
                while (curr && !curr->parent && forests[level].get_degree(curr) == 2 && next &&
                        forests[level].get_degree(next) < 3 && !next->contracts()) {
                    if (!forests[level].try_set_partner_atomic(curr, next)) break;
                    if (forests[level].get_degree(next) == 1) { // If next deg 1 they can combine
                        next->partner = curr;
                        break;
                    }
                    if (!forests[level].try_set_partner_atomic(next, curr)) { // If the CAS fails next was combined from the other side
                        if (next->partner != curr) curr->partner = nullptr;
                        break;
                    }
                    if (next->parent) break; // Stop traversing at a non-root cluster
                    // Get the next two clusters in the chain
                    curr = next->get_other_neighbor(curr);
                    if (curr) next = curr->get_other_neighbor(next);
                }
            }
        } else if (forests[level].get_degree(cluster) == 1) {
            for (auto neighbor : cluster->neighbors) {
                // Combine deg 1 root clusters with deg 1 root or non-root clusters
                if (neighbor && forests[level].get_degree(neighbor) == 1) {
                    forests[level].try_set_partner(cluster, neighbor);
                    forests[level].try_set_partner(neighbor, cluster);
                    break;
                }
                // Combine deg 1 root cluster with deg 2 non-root clusters that don't contract
                if (neighbor && neighbor->parent && forests[level].get_degree(neighbor) == 2) {
                    if (neighbor->contracts()) continue;
                    if (!forests[level].try_set_partner_atomic(neighbor, cluster)) continue;
                    forests[level].try_set_partner(cluster, neighbor);
                    break;
                }
                // Combine deg 1 root cluster with deg 3+ clusters always
                if (neighbor && forests[level].get_degree(neighbor) >= 3) {
                    cluster->partner = neighbor;
                    neighbor->partner = cluster;
                    // If it became high degree, any non-contracting deg 1 clusters must be in the first 2 neighbors
                    if (neighbor->parent)
                    for (auto entry : neighbor->neighbors)
                    if (entry && entry->get_degree() == 1 && entry->parent != neighbor->parent) {
                        entry->partner = neighbor;
                        entry->parent = neighbor->parent;
                    }
                    break;
                }
            }
        }
    });
    // Create new parent clusters for each new combination
    parallel_for(0, R.size(), [&](size_t i) {
        auto cluster = R[i];
        if (cluster->parent && cluster->parent->partner == (Cluster*) 0) return;
        // Make new parent clusters for all new combinations
        auto partner = cluster->partner;
        // Deal with possible multi-combinations
        if (partner && partner->get_degree() >= 3) {
            if (partner->parent) cluster->parent = partner->parent;
            return;
        }
        if (partner && forests[level].get_degree(cluster) >= 3 && !cluster->parent) {
            auto parent = type_allocator<Cluster>::create(
            cluster->priority, default_value);
            parent->partner = (Cluster*) 1;
            cluster->parent = parent;
            for (auto neighbor : cluster->neighbors)
            if (neighbor && neighbor->partner == cluster)
                neighbor->parent = cluster->parent;
            for (auto neighbor : cluster->neighbors_set)
            if (neighbor && neighbor->partner == cluster)
                neighbor->parent = cluster->parent;
            return;
        }
        // Deal with regular pairwise combinations
        if (partner) {
            if (partner->parent && partner->parent != cluster->parent) { // partner is non-root cluster
            cluster->parent = partner->parent;
            partner->partner = nullptr;
            } else if (cluster > partner) { // higher address cluster will do the combination
            auto parent = type_allocator<Cluster>::create(
                cluster->priority + partner->priority, default_value);
            parent->partner = (Cluster*) 1;
            partner->parent = parent;
            cluster->parent = parent;
            }
        } else if (!cluster->parent && forests[level].get_degree(cluster) > 0) { // clusters that don't combine get a new parent
            auto parent = type_allocator<Cluster>::create(
            cluster->priority, default_value);
            parent->partner = (Cluster*) 1;
            cluster->parent = parent;
            cluster->partner = cluster;
        }
    });
    // Fill in the neighbor lists of the new clusters
    parallel_for(0, R.size(), [&](size_t i) {
        auto c1 = R[i];
        auto c2 = c1->partner;
        auto parent = c1->parent;
        if (!c2) return;
        // Deal with possible multi-combinations
        if (c1 != c2 && c2->get_degree() >= 3) return;
        if (c1->get_degree() >= 3) {
            // This means the high degree cluster's parent was newly created
            if (parent->partner == (Cluster*) 1) {
            for (auto neighbor : c1->neighbors)
                if (neighbor && neighbor->parent != parent) {
                parent->insert_neighbor(neighbor->parent);
                if (!neighbor->parent->partner) neighbor->parent->insert_neighbor(parent);
                }
            for (auto neighbor : c1->neighbors_set)
                if (neighbor->parent != parent) {
                parent->insert_neighbor(neighbor->parent);
                if (!neighbor->parent->partner) neighbor->parent->insert_neighbor(parent);
                }
            }
            return;
        }
        // Deal with regular pairwise combinations
        bool new_parent = (parent->partner == (Cluster*) 1);
        if (!new_parent) return;
        if (c1 < c2) return;
        for (int i = 0; i < 3; ++i) {
            if (c1->neighbors[i] && c1->neighbors[i] != c2) {   // Don't add c2's parent (itself)
            parent->insert_neighbor(c1->neighbors[i]->parent);
            if (c1->neighbors[i]->parent->partner == (Cluster*) 0) // Non-root neighbor
                c1->neighbors[i]->parent->insert_neighbor(parent);
            }
        }
        if (c1 != c2)
        for (int i = 0; i < 3; ++i) {
            if (c2->neighbors[i] && c2->neighbors[i] != c1) {   // Don't add c1's parent (itself)
            parent->insert_neighbor(c2->neighbors[i]->parent);
            if (c2->neighbors[i]->parent->partner == (Cluster*) 0) // Non-root neighbor
                c2->neighbors[i]->parent->insert_neighbor(parent);
            }
        }
    });
}

template <typename aug_t>
bool ParallelUFOTree<aug_t>::connected(vertex_t u, vertex_t v) {
    vertex_t root_u = u;
    int level = 0;
    while (forests[level].get_parent(root_u) != -1)
        root_u = forests[level++].get_parent(root_u);
    vertex_t root_v = v;
    level = 0;
    while (forests[level].get_parent(root_v) != -1)
        root_v = forests[level++].get_parent(root_v);
    return root_u == root_v;
}

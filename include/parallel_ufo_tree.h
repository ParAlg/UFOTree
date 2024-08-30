#pragma once

#include "types.h"
#include "util.h"
#include "bridge.h"
#include "hash_bag.h"
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
    // Store R and D, root and del clusters
    hashbag<vertex_t> R;
    hashbag<vertex_t> D;
    hashbag<vertex_t> next_R;
    hashbag<vertex_t> next_D;
    // Helper functions
    void update_tree(sequence<Edge>& updates, bool deletion);
    void recluster_level(int level, bool deletion, sequence<Edge>& U);
    void set_partners(int level);
    void add_parents(int level);
    void add_adjacency(int level);
};

template <typename aug_t>
ParallelUFOTree<aug_t>::ParallelUFOTree(vertex_t n, vertex_t k, QueryType q,
std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d)
    : n(n), query_type(q), f(f), identity(id), default_value(d), R(n), D(n), next_R(n), next_D(n) {
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
    parallel_for(0, updates.size(), [&](size_t i) {
        R.insert(updates[i].src);
        R.insert(updates[i].dst);
        D.insert(forests[0].get_parent(updates[i].src));
        D.insert(forests[0].get_parent(updates[i].dst));
    });
    // Remove parents of low degree level 0 root clusters
    sequence<vertex_t> low_deg = filter(R.extract_all(), [&] (auto v) { return forests[0].get_degree(v) < 2; });
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
    recluster_level(0, deletion, U);
    forests.pop_back();
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_level(int level, bool deletion, sequence<Edge>& U) {
    bool next_R_empty = true;
    bool next_D_empty = true;
    // Delete the edges from our initial updates that are still in level i+1
    if (forests.size() < level) forests.emplace_back();
    if (deletion) forests[level+1].delete_edges(U);

    // Determine the set of clusters that will be deleted in level i+1
    // Collect the new level i root clusters that will be formed by deleting some level i+1 clusters
    R.for_all([&](vertex_t v) {
        vertex_t parent = forests[level].get_parent(v);
        if (parent != NONE) {
            bool low_degree = forests[level+1].get_degree(parent) < 3;
            bool low_fanout = forests[level+1].get_child_count(parent) < 3;
            if (low_degree && low_fanout) {
                auto iter = forests[level].get_neighbor_iterator(v);
                while (vertex_t neighbor = iter->next() != NONE) {
                    if (forests[level].get_parent(neighbor) == parent) {
                        R.insert(neighbor);
                        forests[level].unset_parent(neighbor);
                    }
                }
            }
        }
    });

    // Remove the parents of low degree root clusters
    sequence<vertex_t> low_deg = filter(R.extract_all(), [&] (auto v) { return forests[level].get_degree(v) < 2; });
    sequence<vertex_t> low_deg_parents = parlay::map(low_deg, [&] (auto v) { return forests[level].get_parent(v); });
    forests[level].unset_parents(low_deg);
    forests[level+1].subtract_children(low_deg_parents);

    // Determine the new combinations of the root clusters at this level
    set_partners(level);

    // Get the next set of del clusters and do the deletion of clusters
    D.for_all([&](vertex_t v) {
        vertex_t parent = forests[level+1].get_parent(v);
        if (parent != NONE) {
            next_D.insert(parent);
            next_D_empty = false;
        }
    });
    sequence<vertex_t> del = filter(D.extract_all(), [&] (auto v) {
        bool low_degree = forests[level+1].get_degree(v) < 3;
        bool low_fanout = forests[level+1].get_child_count(v) < 3;
        return (low_degree && low_fanout);
    });
    forests[level+1].delete_vertices(del);

    // Add the new parents and edges to the next level for the set of combinations determined
    add_parents(level);
    add_adjacency(level);

    // Prepare the inputs to the next level
    R.for_all([&](vertex_t v) {
        vertex_t parent = forests[level].get_parent(v);
        if (parent != NONE) {
            next_R.insert(parent);
            next_R_empty = false;
            vertex_t p_parent = forests[level+1].get_parent(parent);
            if (p_parent != NONE) {
                next_D.insert(p_parent);
                next_D_empty = false;
            }
        }
    });
    D.for_all([&](vertex_t v) {
        vertex_t parent = forests[level+1].get_parent(v);
        if (parent != NONE) {
            next_D.insert(parent);
            next_D_empty = false;
        }
    });
    sequence<Edge> next_U = forests[level+1].filter_edges(U);
    // Remove the parents of low degree clusters in the next root clusters
    low_deg = filter(next_R.extract_all(), [&] (auto v) { return forests[level+1].get_degree(v) < 2; });
    low_deg_parents = parlay::map(low_deg, [&] (auto v) { return forests[level+1].get_parent(v); });
    forests[level+1].unset_parents(low_deg);
    if (forests.size() > level+2) forests[level+2].subtract_children(low_deg_parents);
    // Check if the update is complete or recursively recluster the next level
    if (next_R_empty && next_D_empty) return;
    std::swap(R, next_R);
    std::swap(D, next_D);
    next_R.clear();
    next_D.clear();
    recluster_level(level+1, deletion, next_U);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::set_partners(int level) {
    // Perform clustering of the root clusters
    R.for_all([&](vertex_t cluster) {
        if (forests[level].get_parent(cluster)) return;
        // Combine deg 3+ root clusters with deg 1 root or non-root clusters
        if (forests[level].get_degree(cluster) >= 3) {
            // If it became high degree, find any non-contracting deg 1 neighbors
            auto iter = forests[level].get_neighbor_iterator(cluster);
            while (vertex_t neighbor = iter->next() != NONE) {
                if (forests[level].get_degree(neighbor) == 1) {
                    forests[level].try_set_partner(cluster, neighbor);
                    forests[level].try_set_partner(neighbor, cluster);
                    vertex_t parent = forests[level].get_parent(neighbor);
                    if (parent != NONE) {
                        forests[level].unset_parent(neighbor);
                        D.insert(parent);
                    }
                }
            }
        } else if (forests[level].get_degree(cluster) == 2) {
            // Only local maxima in priority with respect to deg 2 clusters will act
            if (!forests[level].is_local_max_priority(cluster)) return;
            // Travel left/right and pair clusters until a deg 3, deg 1, non-root,
            // or combined cluster found
            vertex_t first_neighbor = forests[level].get_first_neighbor(cluster);
            for (auto direction : {0, 1}) {
                auto curr = cluster;
                auto next = direction ? first_neighbor : forests[level].get_other_neighbor(cluster, first_neighbor);
                if (forests[level].has_partner(curr)) {
                    curr = next;
                    next = forests[level].get_other_neighbor(curr, cluster);
                }
                while (curr != NONE && forests[level].get_parent(curr) == NONE && forests[level].get_degree(curr) == 2
                && next != NONE && forests[level].get_degree(next) < 3 && !forests[level].contracts(next)) {
                    if (!forests[level].try_set_partner_atomic(curr, next)) break;
                    if (forests[level].get_degree(next) == 1) { // If next deg 1 they can combine
                        forests[level].try_set_partner(next, curr);
                        break;
                    }
                    if (!forests[level].try_set_partner_atomic(next, curr)) { // If the CAS fails next was combined from the other side
                        if (forests[level].get_partner(next) != curr) forests[level].unset_partner(curr);
                        break;
                    }
                    if (forests[level].get_parent(next) != NONE) break; // Stop traversing at a non-root cluster
                    // Get the next two clusters in the chain
                    curr = forests[level].get_other_neighbor(next, curr);
                    if (curr != NONE) next = forests[level].get_other_neighbor(curr, next);
                }
            }
        } else if (forests[level].get_degree(cluster) == 1) {
            auto iter = forests[level].get_neighbor_iterator(cluster);
            while (vertex_t neighbor = iter->next() != NONE) {
                // Combine deg 1 root clusters with deg 1 root or non-root clusters
                if (neighbor && forests[level].get_degree(neighbor) == 1) {
                    forests[level].try_set_partner(cluster, neighbor);
                    forests[level].try_set_partner(neighbor, cluster);
                    break;
                }
                // Combine deg 1 root cluster with deg 2 non-root clusters that don't contract
                if (forests[level].get_parent(neighbor) != NONE && forests[level].get_degree(neighbor) == 2) {
                    if (forests[level].contracts(neighbor)) continue;
                    if (!forests[level].try_set_partner_atomic(neighbor, cluster)) continue;
                    forests[level].try_set_partner(cluster, neighbor);
                    break;
                }
                // Combine deg 1 root cluster with deg 3+ clusters always
                if (neighbor && forests[level].get_degree(neighbor) >= 3) {
                    forests[level].try_set_partner(cluster, neighbor);
                    forests[level].try_set_partner(neighbor, cluster);
                    // QUINTEN: if a high degree cluster is a non-root cluster, it must have been high degree previously ?
                    // If it became high degree, any non-contracting deg 1 clusters must be in the first 2 neighbors
                    // if (forests[level].get_parent(neighbor) != NONE) {
                    //     auto n_iter = forests[level].get_neighbor_iterator(neighbor);
                    //     while (vertex_t entry = n_iter->next() != NONE) {
                    //         if (forests[level].get_degree(entry) == 1 && forests[level].get_parent(entry) != forests[level].get_parent(neighbor)) {
                    //             forests[level].try_set_partner(entry, neighbor);
                    //             forests[level].set_parent(entry, forests[level].get_parent(neighbor));
                    //         }
                    //     }
                    // }
                    break;
                }
            }
        }
    });
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::add_parents(int level) {
    // Create new parent clusters for each new combination
    R.for_all([&](vertex_t cluster) {
        vertex_t partner = forests[level].get_partner(cluster);
        vertex_t parent = forests[level].get_parent(cluster);
        // Deal with possible multi-combinations
        if (partner != NONE) {
            vertex_t p_partner = forests[level].get_partner(partner);
            vertex_t p_parent = forests[level].get_parent(partner);
            if (forests[level].get_degree(partner) >= 3) {
                if (p_parent != NONE) forests[level].set_parent(cluster, p_parent);
                return;
            }
            if (forests[level].get_degree(cluster) >= 3 && parent == NONE) {
                auto new_parent = cluster;
                forests[level].set_parent(cluster, new_parent);
                auto iter = forests[level].get_neighbor_iterator(cluster);
                while(vertex_t neighbor = iter->next() != NONE)
                    if (forests[level].get_partner(neighbor) == cluster)
                        forests[level].set_parent(neighbor, new_parent);
                return;
            }
        }
        // Deal with regular pairwise combinations
        if (partner != NONE) {
            vertex_t p_partner = forests[level].get_partner(partner);
            vertex_t p_parent = forests[level].get_parent(partner);
            if (p_parent != NONE && p_parent != parent) {   // partner is non-root cluster
                forests[level].set_parent(cluster, p_parent);
                forests[level].unset_partner(partner);
        } else {                                            // two root clusters
                vertex_t new_parent = std::max(cluster, partner);
                forests[level].set_parent(cluster, new_parent);
                forests[level].set_parent(partner, new_parent);
            }
        } else if (parent == NONE && forests[level].get_degree(cluster) > 0) { // clusters that don't combine get a new parent
            vertex_t new_parent = cluster;
            forests[level].set_parent(cluster, new_parent);
            forests[level].try_set_partner(cluster, cluster);
        }
    });
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::add_adjacency(int level) {
    // // Fill in the neighbor lists of the new clusters
    // parallel_for(0, R.size(), [&](size_t i) {
    //     auto c1 = R[i];
    //     auto c2 = c1->partner;
    //     auto parent = c1->parent;
    //     if (!c2) return;
    //     // Deal with possible multi-combinations
    //     if (c1 != c2 && c2->get_degree() >= 3) return;
    //     if (c1->get_degree() >= 3) {
    //         // This means the high degree cluster's parent was newly created
    //         if (parent->partner == (Cluster*) 1) {
    //         for (auto neighbor : c1->neighbors)
    //             if (neighbor && neighbor->parent != parent) {
    //             parent->insert_neighbor(neighbor->parent);
    //             if (!neighbor->parent->partner) neighbor->parent->insert_neighbor(parent);
    //             }
    //         for (auto neighbor : c1->neighbors_set)
    //             if (neighbor->parent != parent) {
    //             parent->insert_neighbor(neighbor->parent);
    //             if (!neighbor->parent->partner) neighbor->parent->insert_neighbor(parent);
    //             }
    //         }
    //         return;
    //     }
    //     // Deal with regular pairwise combinations
    //     bool new_parent = (parent->partner == (Cluster*) 1);
    //     if (!new_parent) return;
    //     if (c1 < c2) return;
    //     for (int i = 0; i < 3; ++i) {
    //         if (c1->neighbors[i] && c1->neighbors[i] != c2) {   // Don't add c2's parent (itself)
    //         parent->insert_neighbor(c1->neighbors[i]->parent);
    //         if (c1->neighbors[i]->parent->partner == (Cluster*) 0) // Non-root neighbor
    //             c1->neighbors[i]->parent->insert_neighbor(parent);
    //         }
    //     }
    //     if (c1 != c2)
    //     for (int i = 0; i < 3; ++i) {
    //         if (c2->neighbors[i] && c2->neighbors[i] != c1) {   // Don't add c1's parent (itself)
    //         parent->insert_neighbor(c2->neighbors[i]->parent);
    //         if (c2->neighbors[i]->parent->partner == (Cluster*) 0) // Non-root neighbor
    //             c2->neighbors[i]->parent->insert_neighbor(parent);
    //         }
    //     }
    // });
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

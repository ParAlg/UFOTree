#pragma once
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include "types.h"
#include "util.h"
#include "bridge.h"
#include "hash_bag.h"
#include "parallel_ufo_cluster.h"


namespace dgbs {

template <typename aug_t = empty_t>
class ParallelUFOTree {
    using Cluster = ParallelUFOCluster<aug_t>;
    using allocator = parlay::type_allocator<Cluster>;
public:
    // UFO tree interface
    ParallelUFOTree(vertex_t n, vertex_t k);
    ~ParallelUFOTree();
    void batch_link(parlay::sequence<std::pair<int, int>>& links);
    void batch_cut(parlay::sequence<std::pair<int, int>>& cuts);
    bool connected(vertex_t u, vertex_t v);
    // Testing helpers
    bool is_valid();
    int get_height(vertex_t v);
    void print_tree();
private:
    // Class data and parameters
    std::vector<Cluster> leaves;
    // Helper functions
    void recluster_tree(parlay::sequence<std::pair<int, int>>& updates, UpdateType update_type);
};

template <typename aug_t>
ParallelUFOTree<aug_t>::ParallelUFOTree(vertex_t n, vertex_t k) : leaves(n) {}

template <typename aug_t>
ParallelUFOTree<aug_t>::~ParallelUFOTree() {
    allocator::finish();
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_link(parlay::sequence<std::pair<int, int>>& links) {
    recluster_tree(links, INSERT);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_cut(parlay::sequence<std::pair<int, int>>& cuts) {
    recluster_tree(cuts, DELETE);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_tree(parlay::sequence<std::pair<int, int>>& updates, UpdateType update_type) {
    // ==============
    // INITIALIZATION
    // ==============
    parlay::sequence<Cluster*> root_clusters;
    parlay::sequence<Cluster*> del_clusters;
    parlay::sequence<std::pair<Cluster*, Cluster*>> dir_updates;

    // The intial dir updates are just the batch of links or cuts in both directions.
    dir_updates = parlay::tabulate(2*updates.size(), [&] (size_t i) {
        if (i % 2 == 0) return std::make_pair(&leaves[updates[i].first], &leaves[updates[i].second]);
        return std::make_pair(&leaves[updates[i].second], &leaves[updates[i].first]);
    });
    
    // Group by endpoint. Then group by parent.
    auto dir_update_groups = parlay::group_by_key(dir_updates);
    auto parents = parlay::delayed_tabulate(dir_update_groups.size(), [&] (size_t i) {
        return std::make_pair(dir_update_groups[i].first->parent, dir_update_groups[i].first);
    });
    auto parent_groups = parlay::group_by_key(parents);

    // The initial root clusters are just all the vertices that got an update.
    root_clusters = parlay::map(dir_update_groups, [&] (auto group) {
        return group.first;
    });

    // The intial del clusters are just all parents of vertices that got an update.
    del_clusters = parlay::map_maybe(parent_groups, [&] (auto group) -> std::optional<Cluster*> {
        if (!group.first) return std::nullopt;
        return group.first;
    });

    /* This is the main loop over all levels of the tree. At each level i, we assume that `root_clusters`is
    populated with the level i clusters formed from contractions in the previous level, and that `del_clusters`
    is populated with the level i+1 clusters that were parents of del clusters at the previous level, along
    with the parents of any nodes that started a new remove ancestor path due to a root cluster merging into a
    non-root cluster that didn't already contract.
    (1) First we take the level i+1 `del_clusters`, and populate `next_del_clusters` by mapping them to their
    level i+2 parents. We also determine which ones will actually be deleted. During this step we delete the
    del clusters. We also find any level i+1 non-del children of the level i+2 clusters which will be deleted.
    These will become `next_root_clusters` at the next level i+1.
    (2) For del clusters that don't get deleted, we need to update their adjacency lists based on the batch of
    links or cuts. For cuts it is necessary to do this after determining the `next_root_clusters` in step 1, so
    that the siblings of a cluster can still be found.
    (3) Next we, will determine the new set of contractions over the level i root clusters using an MIS. While
    doing this we also allocate new level i+1 parent nodes and update the parent pointers of the root clusters.
    These new parents become the `next_root_clusters` along with those determined in step 1. Any contraction
    with a non-root cluster adds its level i+2 grandparent to `next_del_clusters`. After, we fill the adjacency
    lists of our new parent clusters and add them to the adjacency list of their neighbors. */
    while (root_clusters.size() > 0 || del_clusters.size() > 0) {
        // =======
        // PHASE 1
        // =======

        // Group del clusters by parent.
        auto parents = parlay::delayed_tabulate(del_clusters.size(), [&] (size_t i) {
            return std::make_pair(del_clusters[i]->parent, del_clusters[i]);
        });
        auto parent_groups = parlay::group_by_key(parents);
        // Get the next del clusters.
        auto next_del_clusters = parlay::delayed::map_maybe(parent_groups, [&] (auto group) -> std::optional<Cluster*> {
            if (!group.first) return std::nullopt;
            return group.first;
        });
        // Get next root clusters, clear their parent pointers, delete some current del clusters.
        auto next_root_clusters = parlay::delayed::flatten(parlay::delayed_tabulate(parent_groups.size(), [&] (size_t i) {
            auto& [parent, children] = parent_groups[i];
            parlay::sequence<Cluster*> local_root_clusters;
            if (parent == nullptr) return local_root_clusters;

            Cluster* max = *parlay::max_element(children, [&] (Cluster* x, Cluster* y) { return x->get_degree() < y->get_degree(); });
            size_t max_degree = max->get_degree();

            if (max_degree == 1) {
                local_root_clusters = std::move(children);
                Cluster* center = max->get_neighbor();
                size_t center_degree = center->get_degree();
                if (center_degree - children.size() < 3) {
                    center->delete_neighbors(children); // Should replace delete / reinsert with parallel iteration
                    Cluster* neighbor1 = center->get_neighbor();
                    if (neighbor1->parent == parent) local_root_clusters.push_back(neighbor1);
                    if (center_degree - children.size() == 2) {
                        Cluster* neighbor2 = center->get_other_neighbor(neighbor1);
                        if (neighbor2->parent == parent) local_root_clusters.push_back(center->get_other_neighbor(neighbor2));
                    }
                    center->insert_neighbors(children);
                    allocator::destroy(parent);
                }
            }

            else if (max_degree == 2) {
                local_root_clusters.push_back(max);
                Cluster* neighbor = max->get_neighbor();
                if (neighbor->parent = parent) local_root_clusters.push_back(neighbor);
                else local_root_clusters.push_back(max->get_other_neighbor(neighbor));
                allocator::destroy(parent);
            }

            else if (max_degree >= 3) {
                local_root_clusters = std::move(children);
                if (max_degree - (children.size()-1) < 3) {
                    max->delete_neighbors(children); // Should replace delete / reinsert with parallel iteration
                    Cluster* neighbor1 = max->get_neighbor();
                    if (neighbor1->parent == parent) local_root_clusters.push_back(neighbor1);
                    if (max_degree - children.size() == 2) {
                        Cluster* neighbor2 = max->get_other_neighbor(neighbor1);
                        if (neighbor2->parent == parent) local_root_clusters.push_back(max->get_other_neighbor(neighbor2));
                    }
                    max->insert_neighbors(children);
                    allocator::destroy(parent);
                } else {
                    local_root_clusters.erase(parlay::find(local_root_clusters, max));
                }
            }

            parlay::parallel_for(0, local_root_clusters.size(), [&] (size_t i) {
                local_root_clusters[i]->parent = nullptr;
            });
            return local_root_clusters;
        }));
        // =======
        // PHASE 2
        // =======

        // Insert or delete the edges at this level.
        auto dir_update_groups = parlay::group_by_key(dir_updates);
        parlay::parallel_for(0, dir_update_groups.size(), [&] (size_t i) {
            auto& [cluster, neighbors] = dir_update_groups[i];
            if (update_type == INSERT) cluster->insert_neighbors(neighbors);
            else cluster->delete_neighbors(neighbors);
        });

        // Map the dir updates to the next level.
        dir_updates = parlay::map_maybe(dir_updates, [&] (std::pair<Cluster*, Cluster*> x) -> std::optional<std::pair<Cluster*, Cluster*>> {
            if (x.first->parent != x.second->parent)
                return std::make_pair(x.first->parent, x.second->parent);
            return std::nullopt;
        });

        // =======
        // PHASE 3
        // =======

        /* We also need to handle the case where a root cluster became high degree because of many
        insertions, but it needs to find a degree 1 neighbor that is a non-root cluster because it
        previously happened to not contract with this cluster. Fortunately we can detect when this
        will happen by looking at the number of edges inserted, and we can grab that deg 1 neighbor
        easily before inserting all the edges and making this cluster high degree.
        This might mean populating the adjacency lists should happen during the phase 3 of the next
        level up. */


        // We can probably eliminate the extra `partner` field per cluster by using the `parent` field smartly.
        parlay::parallel_for(0, root_clusters.size(), [&] (size_t i) {
            Cluster* cluster = root_clusters[i];
            if (cluster->get_degree() == 1) {
                Cluster* neighbor = cluster->get_neighbor();
                // Combine deg 1 root clusters with deg 1 root or non-root clusters
                if (neighbor->get_degree() == 1) {
                    cluster->partner = neighbor;
                    neighbor->partner = cluster;
                }
                // Combine deg 1 root cluster with deg 2 non-root clusters that don't contract
                else if (neighbor->get_degree() == 2 && neighbor->parent) {
                    if (!neighbor->contracts())
                        if (!neighbor->partner && gbbs::CAS(&neighbor->partner, (Cluster*) nullptr, cluster))
                            cluster->partner = neighbor;
                }
                // Combine deg 1 root cluster with deg 3+ clusters always
                else if (neighbor->get_degree() >= 3) {
                    cluster->partner = neighbor;
                    neighbor->partner = cluster;
                }
            }
            else if (cluster->get_degree() == 2) {
                // Only local maxima in priority with respect to deg 2 clusters will act
                Cluster* neighbor1 = cluster->get_neighbor();
                Cluster* neighbor2 = cluster->get_other_neighbor(neighbor1);
                uint64_t hash = hash64((uintptr_t) cluster);
                uint64_t hash1 = hash64((uintptr_t) neighbor1);
                uint64_t hash2 = hash64((uintptr_t) neighbor2);
                if (neighbor1->get_degree() == 2 && (hash1 > hash || (hash1 == hash && neighbor1 > cluster))) return;
                if (neighbor2->get_degree() == 2 && (hash2 > hash || (hash2 == hash && neighbor2 > cluster))) return;
                // Travel left/right and pair clusters until a deg 3, deg 1, non-root, or partnered cluster is found
                for (auto direction : {0, 1}) {
                    auto curr = cluster;
                    auto next = direction ? neighbor1 : neighbor2;
                    if (curr->partner) {
                        curr = next;
                        next = curr->get_other_neighbor(cluster);
                    }
                    while (curr && !curr->parent && curr->get_degree() == 2 && next && next->get_degree() < 3 && !next->contracts()) {
                        if (curr->partner || !gbbs::CAS(&curr->partner, (Cluster*) nullptr, next)) return;
                        if (next->get_degree() == 1) { // If next deg 1 they can combine
                            next->partner = curr;
                            return;
                        }
                        if (next->partner || !gbbs::CAS(&next->partner, (Cluster*) nullptr, curr)) { // If the CAS fails next was combined from the other side
                            if (next->partner != curr) curr->partner = nullptr;
                            return;
                        }
                        if (next->parent) return; // Stop traversing at a non-root cluster
                        // Get the next two clusters in the chain
                        curr = next->get_other_neighbor(curr);
                        if (curr) next = curr->get_other_neighbor(next);
                    }
                }
            }
        });
    }
}

template <typename aug_t>
bool ParallelUFOTree<aug_t>::connected(vertex_t u, vertex_t v) {
    return (leaves[u].get_root() == leaves[v].get_root());
}

}

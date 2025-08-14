#pragma once
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include "types.h"
#include "util.h"
#include "bridge.h"
#include "hash_bag.h"
#include "parallel_ufo_cluster.h"


extern parlay::internal::timer timer1;
extern parlay::internal::timer timer2;
extern parlay::internal::timer timer3;
extern parlay::internal::timer timer4;
extern parlay::internal::timer timer5;

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
// private:
    // Class data and parameters
    std::vector<Cluster> leaves;
    // Helper functions
    void recluster_tree(parlay::sequence<std::pair<int, int>>& updates, UpdateType update_type);

    static parlay::sequence<Cluster*> recluster_root_clusters(parlay::sequence<Cluster*>& root_clusters);
    static Cluster* recluster_degree_one_root(Cluster* root_cluster);
    static parlay::sequence<Cluster*> recluster_degree_two_root(Cluster* root_cluster);
    static Cluster* recluster_high_degree_root(Cluster* root_cluster);
    static parlay::sequence<Cluster*> create_new_parents(parlay::sequence<Cluster*>& root_clusters);
    static inline bool is_local_max(Cluster* c);
};

template <typename aug_t>
ParallelUFOTree<aug_t>::ParallelUFOTree(vertex_t n, vertex_t k) : leaves(n) {}

template <typename aug_t>
ParallelUFOTree<aug_t>::~ParallelUFOTree() {
    auto clusters_to_delete = parlay::flatten(parlay::tabulate(leaves.size(), [&] (size_t i) {
        parlay::sequence<Cluster*> clusters;
        Cluster* curr = leaves[i].parent;
        while (curr && curr != (Cluster*) 1) {
            Cluster* next = AtomicLoad(&curr->parent);
            if (next != (Cluster*) 1 && CAS(&curr->parent, next, (Cluster*) 1))
                clusters.push_back(curr);
            else break;
            curr = next;
        }
        return clusters;
    }));
    parlay::parallel_for(0, clusters_to_delete.size(), [&] (size_t i) {
        allocator::destroy(clusters_to_delete[i]);
    });
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
        if (i % 2 == 0) return std::make_pair(&leaves[updates[i/2].first], &leaves[updates[i/2].second]);
        return std::make_pair(&leaves[updates[i/2].second], &leaves[updates[i/2].first]);
    });

    // Group by endpoint. Then group by parent.
    auto dir_update_groups = integer_group_by_key_inplace(dir_updates);
    auto parents = parlay::map(dir_update_groups, [&] (auto group) {
        return std::make_pair(group.first->parent, group.first);
    });
    auto parent_groups = integer_group_by_key_inplace(parents);

    // The level 0 root clusters are children of level 1 parents that will get deleted, or deg 1 endpoints.
    auto unflattened_root_clusters = parlay::tabulate(parent_groups.size(), [&] (size_t i) {
        auto& [parent, children] = parent_groups[i];
        if (parent == nullptr) return parlay::map(children, [&] (auto child) { return child.second; });
        parlay::sequence<Cluster*> local_root_clusters;

        Cluster* max = (*parlay::max_element(children, [&] (auto x, auto y) { return x.second->get_degree() < y.second->get_degree(); })).second;
        size_t max_degree = max->get_degree();

        if (max_degree == 1) {
            Cluster* center = max->get_neighbor();
            if (center->parent != max->parent) {
                local_root_clusters.push_back(max);
                parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
            } else if (center->get_degree() < children.size() + 3) {
                local_root_clusters = center->filter_neighbors([&] (auto neighbor) {
                    return neighbor->parent == parent;
                });
                local_root_clusters.push_back(center);
                parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
            }
        }

        else {
            if (max_degree < (children.size()-1) + 3) {
                local_root_clusters = max->filter_neighbors([&] (auto neighbor) {
                    return neighbor->parent == parent;
                });
                local_root_clusters.push_back(max);
                parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
            } else {
                local_root_clusters = parlay::map_maybe(children, [&] (auto child) -> std::optional<Cluster*> {
                    if (child.second != max) return child.second;
                    return std::nullopt;
                });
            }
        }

        return local_root_clusters;
    });
    root_clusters = parlay::flatten(unflattened_root_clusters);

    // The intial del clusters are just all parents of vertices that got an update.
    del_clusters = parlay::filter(parlay::delayed_tabulate(parent_groups.size(), [&] (size_t i) {
        // Null the parents of all root clusters
        if (unflattened_root_clusters[i].size() < 100) {
            for (int j = 0; j < unflattened_root_clusters[i].size(); j++) {
                unflattened_root_clusters[i][j]->parent = nullptr;
            }
        } else {
            parlay::parallel_for(0, unflattened_root_clusters[i].size(), [&] (size_t j) {
                unflattened_root_clusters[i][j]->parent = nullptr;
            });
        }
        // Return the parent
        return parent_groups[i].first;
    }), [&] (auto cluster) { return cluster != nullptr; });


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
    size_t level = 0;
    while (root_clusters.size() > 0 || del_clusters.size() > 0) {
        // =======
        // PHASE 2
        // =======

        // Insert or delete the edges at this level.
        parlay::parallel_for(0, dir_update_groups.size(), [&] (size_t i) {
            auto& [cluster, edges] = dir_update_groups[i];
            auto neighbors = parlay::map(edges, [&] (auto x) { return x.second; });
            if (update_type == INSERT) cluster->insert_neighbors_sorted(neighbors);
            else cluster->delete_neighbors_sorted(neighbors);
        });

        // Map the dir updates to the next level.
        dir_updates = parlay::map_maybe(dir_updates, [&] (std::pair<Cluster*, Cluster*> x) -> std::optional<std::pair<Cluster*, Cluster*>> {
            if (x.first->parent && x.second->parent && x.first->parent != x.second->parent)
                return std::make_pair(x.first->parent, x.second->parent);
            return std::nullopt;
        });
        dir_update_groups = integer_group_by_key_inplace(dir_updates);

        // =======
        // PHASE 3
        // =======

        /* We also need to handle the case where a root cluster became high degree because of many
        insertions, but it needs to find a degree 1 neighbor that is a non-root cluster because it
        previously happened to not contract with this cluster. Fortunately we can detect when this
        will happen by looking at the number of edges inserted, and we can grab that deg 1 neighbor
        easily before inserting all the edges and making this cluster high degree.
        This means populating the adjacency lists should happen during phase 3 of the next level. */

        auto additional_del_clusters = recluster_root_clusters(root_clusters);

        // =======
        // PHASE 1
        // =======

        del_clusters.append(additional_del_clusters);

        // Group del clusters by parent.
        auto parents = parlay::tabulate(del_clusters.size(), [&] (size_t i) {
            return std::make_pair(del_clusters[i]->parent, del_clusters[i]);
        });
        auto parent_groups = integer_group_by_key_inplace(parents);
        // Get the next del clusters.
        auto next_del_clusters = parlay::map_maybe(parent_groups, [&] (auto group) -> std::optional<Cluster*> {
            if (!group.first) return std::nullopt;
            return group.first;
        });
        // Get next root clusters, clear their parent pointers, delete some current del clusters.
        // This should return any del cluster that doesn't get deleted, but whose parent get's deleted.
        // This should also return any cluster that is not a del cluster, but is a child of a deleted parent.
        auto next_root_clusters2 = parlay::flatten(parlay::tabulate(parent_groups.size(), [&] (size_t i) {
            auto& [parent, children] = parent_groups[i];
            parlay::sequence<Cluster*> local_root_clusters;
            if (parent == nullptr) return local_root_clusters;

            Cluster* max = (*parlay::max_element(children, [&] (auto x, auto y) { return x.second->get_degree() < y.second->get_degree(); })).second;
            size_t max_degree = max->get_degree();

            if (max_degree == 1) {
                Cluster* center = max->get_neighbor();
                if (center->parent != max->parent) {
                    if (!max->partner) local_root_clusters.push_back(max);
                    parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
                } else if (center->get_degree() < children.size() + 3) {
                    local_root_clusters = center->filter_neighbors([&] (auto neighbor) {
                        return !neighbor->partner && neighbor->parent == parent;
                    });
                    if (!center->partner) local_root_clusters.push_back(center);
                    parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
                } else {
                    local_root_clusters = parlay::map_maybe(children, [&] (auto x) -> std::optional<Cluster*> {
                        if (!x.second->partner) return x.second;
                        return std::nullopt;
                    });
                }
            }

            else if (max_degree == 2) {
                if (!max->partner) local_root_clusters.push_back(max);
                Cluster* neighbor1 = max->get_neighbor();
                Cluster* neighbor2 = max->get_other_neighbor(neighbor1);
                if (!neighbor1->partner && neighbor1->parent == parent) local_root_clusters.push_back(neighbor1);
                else if (!neighbor2->partner && neighbor2->parent == parent) local_root_clusters.push_back(neighbor2);
                parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
            }

            else if (max_degree >= 3) {
                if (max_degree < (children.size()-1) + 3) {
                    local_root_clusters = max->filter_neighbors([&] (auto neighbor) {
                        return !neighbor->partner && neighbor->parent == parent;
                    });
                    if (!max->partner) local_root_clusters.push_back(max);
                    parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
                } else {
                    local_root_clusters = parlay::map_maybe(children, [&] (auto x) -> std::optional<Cluster*> {
                        if (!x.second->partner && x.second != max) return x.second;
                        return std::nullopt;
                    });
                }
            }

            return local_root_clusters;
        }));

        // Clear pointers to nodes that will get deleted.
        auto del_edges = parlay::flatten(parlay::tabulate(parent_groups.size(), [&] (size_t i) {
            auto& [parent, children] = parent_groups[i];
            auto local_edges = parlay::flatten(parlay::tabulate(children.size(), [&] (size_t j) {
                parlay::sequence<std::pair<Cluster*, Cluster*>> local_local_edges;
                if (children[j].second->partner) {
                    local_local_edges = parlay::map(children[j].second->filter_neighbors([&] (auto neighbor) {
                        return !neighbor->partner;
                    }), [&] (auto cluster) {
                        return std::make_pair(cluster, children[j].second);
                    });
                }
                return local_local_edges;
            }));
            return local_edges;
        }));
        auto del_edge_groups = integer_group_by_key_inplace(del_edges);

        parlay::parallel_for(0, del_edge_groups.size(), [&] (size_t i) {
            auto& [cluster, edges] = del_edge_groups[i];
            auto neighbors = parlay::map(edges, [&] (auto x) { return x.second; });
            cluster->delete_neighbors(neighbors);
        });

        // Delete clusters.
        parlay::parallel_for(0, del_clusters.size(), [&] (size_t i) {
            if (del_clusters[i]->partner) allocator::destroy(del_clusters[i]);
        });

        // =======
        // PHASE 3
        // =======

        // This returns only the new clusters that were created during the reclustering at this level.
        auto next_root_clusters1 = create_new_parents(root_clusters);

        // Fill the adjacency lists of new clusters
        // This should work well for only low degree cases
        parlay::parallel_for(0, root_clusters.size(), [&] (size_t i) {
            auto cluster = root_clusters[i];
            if (!cluster->parent) return; // Only deg 0

            // Clear partner pointers
            Cluster* partner = AtomicLoad(&cluster->partner);
            if (partner) {
                Cluster* partner_partner = AtomicLoad(&partner->partner);
                if (partner_partner != cluster) { // Non-root partner
                    AtomicStore(&partner->partner, (Cluster*) nullptr);
                    AtomicStore(&cluster->partner, (Cluster*) nullptr);
                } else if (cluster < partner) { // Tie-break
                    AtomicStore(&partner->partner, (Cluster*) nullptr);
                    AtomicStore(&cluster->partner, (Cluster*) nullptr);
                }
            }

            // Fill adjacency lists at the next level up
            cluster->for_all_neighbors([&] (auto neighbor) {
                if (neighbor->parent != cluster->parent) {
                    cluster->parent->insert_neighbor(neighbor->parent);
                    neighbor->parent->insert_neighbor(cluster->parent);
                }
            });
        });

        // ===============
        // PREP NEXT LEVEL
        // ===============
        parlay::parallel_for(0, next_root_clusters2.size(), [&] (size_t i) {
            next_root_clusters2[i]->parent = nullptr;
        });
        root_clusters = parlay::append(next_root_clusters1, next_root_clusters2);
        del_clusters = std::move(next_del_clusters);
        level++;
    }
}

template <typename aug_t>
parlay::sequence<ParallelUFOCluster<aug_t>*> ParallelUFOTree<aug_t>::recluster_root_clusters(parlay::sequence<Cluster*>& root_clusters) {
    // This function sets the partner fields for all root clusters. For a non-root
    // cluster, we mark its partner field with 0x1. For root clusters that don't
    // combine with anything, we leave its partner field empty. For high degree
    // root clusters, we don't assign it a partner, but we add a parent for it in
    // this part. All other root clusters receive no parent at this point.
    return parlay::flatten(parlay::tabulate(root_clusters.size(), [&] (size_t i) {
        parlay::sequence<Cluster*> del_clusters;
        Cluster* cluster = root_clusters[i];
        if (cluster->get_degree() == 1) {
            Cluster* del_cluster = recluster_degree_one_root(cluster);
            if (del_cluster) del_clusters.push_back(del_cluster);
        }
        else if (cluster->get_degree() == 2) {
            del_clusters = recluster_degree_two_root(cluster);
        }
        else if (cluster->get_degree() >= 3) {
            Cluster* del_cluster = recluster_high_degree_root(cluster);
            if (del_cluster) del_clusters.push_back(del_cluster);
        }
        return del_clusters;
    }));
}

template <typename aug_t>
ParallelUFOCluster<aug_t>* ParallelUFOTree<aug_t>::recluster_degree_one_root(Cluster* cluster) {
    // Always partner with a degree 1 or 3+ neighbor. For degree 2
    // neighbors, only attempt to partner with it if it is a non-root
    // cluster that does not already contract. Partnering with degree
    // 2 root clusters will be handled from the deg 2 cluster's side.
    // Return the parent of a non-root combination, if any.
    Cluster* neighbor = cluster->get_neighbor();
    if (neighbor->get_degree() == 1) { // Combine deg 1 root clusters with deg 1 root or non-root clusters
        cluster->partner = neighbor;
        if (neighbor->parent) {
            neighbor->partner = (Cluster*) 1; // Mark non-root contracting cluster's partner field
            return neighbor->parent;
        }
    }
    else if (neighbor->get_degree() == 2) { // Combine deg 1 root cluster with deg 2 non-root clusters that don't contract
        if (neighbor->parent && !neighbor->contracts()) {
            if (CAS(&neighbor->partner, (Cluster*) nullptr, (Cluster*) 1)) { // Mark non-root contracting cluster's partner field
                cluster->partner = neighbor;
                return neighbor->parent;
            }
        }
    }
    else { // Combine deg 1 root cluster with possible deg 3+ non-root clusters
        cluster->partner = neighbor;
    }
    return nullptr;
}

template <typename aug_t>
parlay::sequence<ParallelUFOCluster<aug_t>*> ParallelUFOTree<aug_t>::recluster_degree_two_root(Cluster* cluster) {
    parlay::sequence<Cluster*> local_del_clusters;
    // Only local maxima in priority with respect to deg 2 clusters will act
    if (!is_local_max(cluster)) return local_del_clusters;
    // Travel left/right and pair clusters until a deg 3+, deg 1, non-root, or partnered cluster is found
    Cluster* neighbor1 = cluster->get_neighbor();
    Cluster* neighbor2 = cluster->get_other_neighbor(neighbor1);
    for (bool direction : {0, 1}) {
        Cluster* curr = cluster;
        Cluster* next = direction ? neighbor1 : neighbor2;
        if (AtomicLoad(&curr->partner)) {
            curr = next;
            next = curr->get_other_neighbor(cluster);
        }
        while (curr && !curr->parent && curr->get_degree() == 2 && next && next->get_degree() < 3 && !next->contracts()) {
            if (!CAS(&curr->partner, (Cluster*) nullptr, next)) break;
            if (next->get_degree() == 1) { // If next deg 1 they can definitely combine
                if (!next->parent) next->partner = curr;
                else next->partner = (Cluster*) 1; // Mark non-root contracting cluster's partner field
            } else {
                Cluster* new_partner = next->parent ? (Cluster*) 1 : curr; // Mark non-root contracting cluster's partner field
                if (!CAS(&next->partner, (Cluster*) nullptr, new_partner)) { // If the CAS fails next was combined from the other side
                    if (AtomicLoad(&next->partner) != curr) // Other side combined with the opposite cluster (you got left hanging)
                        AtomicStore(&curr->partner, (Cluster*) nullptr);
                    break;
                }
            }
            // Both CAS's succeeded or next was degree 1
            if (next->parent) { // Stop traversing at a non-root cluster
                local_del_clusters.push_back(next->parent); // This happens at most twice, once per direction
                break;
            }
            if (next->get_degree() == 1) // Stop traversing at deg 1 cluster
                break;
            // Get the next two clusters in the chain
            curr = next->get_other_neighbor(curr);
            if (curr) next = curr->get_other_neighbor(next);
        }
    }
    return local_del_clusters;
}

template <typename aug_t>
ParallelUFOCluster<aug_t>* ParallelUFOTree<aug_t>::recluster_high_degree_root(Cluster* cluster) {
    // Create the new parent for a high degree root cluster.
    // Find at most one possible non-root degree 1 neighbor,
    // combine with it, and return its parent as del cluster.
    Cluster* del_cluster = nullptr;
    Cluster* parent = allocator::create();
    cluster->parent = parent;
    cluster->for_all_neighbors([&] (auto neighbor) {
        if (neighbor->get_degree() == 1 && neighbor->parent && neighbor->parent != parent) {
            neighbor->parent->partner = (Cluster*) 1;
            del_cluster = neighbor->parent; // Only once in this loop
            neighbor->parent = parent;
        }
    });
    return del_cluster;
}

template <typename aug_t>
parlay::sequence<ParallelUFOCluster<aug_t>*> ParallelUFOTree<aug_t>::create_new_parents(parlay::sequence<Cluster*>& root_clusters) {
    // Only returns the brand new clusters to be root clusters at the next level,
    // not the parents of non-root clusters. Those may become root clusters also,
    // but it will be determined by the later step which checks if the grandparent
    // should be deleted.
    return parlay::map_maybe(root_clusters, [&] (Cluster* cluster) -> std::optional<Cluster*> {
        if (cluster->get_degree() == 0) return std::nullopt;
        if (cluster->get_degree() >= 3) return cluster->parent;
        Cluster* partner = cluster->partner;
        if (partner) {
            if (partner->partner != cluster) { // Non-root partner
                cluster->parent = partner->parent;
                return std::nullopt;
            } else if (cluster < partner) { // Tie-break for two partnered root clusters
                Cluster* parent = allocator::create();
                cluster->parent = parent;
                partner->parent = parent;
                return parent;
            } else { // Lost tie-break, only return new parent once
                return std::nullopt;
            }
        }
        else { // Non-combining root cluster gets its own parent
            Cluster* parent = allocator::create();
            cluster->parent = parent;
            return parent;
        }
    });
}

template <typename aug_t>
bool ParallelUFOTree<aug_t>::is_local_max(Cluster* c) {
    // Assumes the input is a degree 2 cluster
    Cluster* neighbor1 = c->get_neighbor();
    Cluster* neighbor2 = c->get_other_neighbor(neighbor1);
    uint64_t hash = hash64((uintptr_t) c);
    uint64_t hash1 = hash64((uintptr_t) neighbor1);
    uint64_t hash2 = hash64((uintptr_t) neighbor2);
    if (neighbor1->get_degree() == 2 && (!AtomicLoad(&neighbor1->parent) || AtomicLoad(&neighbor1->partner)))
        if (hash1 > hash || (hash1 == hash && neighbor1 > c))
            return false;
    if (neighbor2->get_degree() == 2 && (!AtomicLoad(&neighbor2->parent) || AtomicLoad(&neighbor2->partner)))
        if (hash2 > hash || (hash2 == hash && neighbor2 > c))
            return false;
    return true;
}

template <typename aug_t>
bool ParallelUFOTree<aug_t>::connected(vertex_t u, vertex_t v) {
    return (leaves[u].get_root() == leaves[v].get_root());
}

}

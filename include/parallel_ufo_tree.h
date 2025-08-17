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
    static constexpr uintptr_t NULL_PTR = 0;
    static constexpr uintptr_t DEL_MARK = 1;
    static constexpr uintptr_t NON_ROOT_MARK = 2;
    static constexpr uintptr_t NEW_PAR_MARK = 3;
    using allocator = parlay::type_allocator<Cluster>;
    using EdgeSlice = parlay::slice<std::pair<Cluster*, Cluster*>*, std::pair<Cluster*, Cluster*>*>;
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

    static parlay::sequence<Cluster*> process_initial_clusters(parlay::sequence<std::pair<Cluster*, EdgeSlice>>& parent_groups);
    static parlay::sequence<Cluster*> process_del_clusters(parlay::sequence<std::pair<Cluster*, EdgeSlice>>& parent_groups);
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
        while (curr && curr != (Cluster*) DEL_MARK) {
            Cluster* next = AtomicLoad(&curr->parent);
            if (next != (Cluster*) DEL_MARK && CAS(&curr->parent, next, (Cluster*) DEL_MARK))
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
    auto dir_update_groups = integer_group_by_key_inplace(dir_updates);

    // Group the affected vertices by parent.
    auto parents = parlay::map(dir_update_groups, [&] (auto group) {
        // For deletion batches, compute new degrees of del clusters.
        if (update_type == DELETE) {
            group.first->degree = group.first->get_degree() - group.second.size();
        }
        return std::make_pair(group.first->parent, group.first);
    });
    auto parent_groups = integer_group_by_key_inplace(parents);

    // For deletion batches, map the level 0 updates to level 1 before processing level 0.
    parlay::sequence<std::pair<Cluster*, Cluster*>> temp_dir_updates;
    if (update_type == DELETE) {
        temp_dir_updates = parlay::map_maybe(dir_updates, [&] (std::pair<Cluster*, Cluster*> x) -> std::optional<std::pair<Cluster*, Cluster*>> {
            if (x.first->parent && x.second->parent && x.first->parent != x.second->parent)
                return std::make_pair(x.first->parent, x.second->parent);
            return std::nullopt;
        });
    }

    // The initial root clusters are children of level 1 parents that will get deleted, or deg 1 endpoints.
    root_clusters = process_initial_clusters(parent_groups);
    parlay::parallel_for(0, root_clusters.size(), [&] (size_t i) {
        root_clusters[i]->parent = nullptr;
    });

    // For deletion batches, do the deletions at level 0, finish mapping the updates to level 1.
    if (update_type == DELETE) {
        parlay::parallel_for(0, dir_update_groups.size(), [&] (size_t i) {
            auto& [cluster, edges] = dir_update_groups[i];
            auto neighbors = parlay::map(edges, [&] (auto x) { return x.second; });
            cluster->delete_neighbors(neighbors);
        });

        dir_updates = std::move(temp_dir_updates);
        dir_update_groups = integer_group_by_key_inplace(dir_updates);
    }

    // The intial del clusters are just all parents of vertices that got an update.
    del_clusters = parlay::map_maybe(parent_groups, [&] (auto group) -> std::optional<Cluster*> {
        if (group.first) return group.first;
        return std::nullopt;
    });


    /* This is the main loop over all levels of the tree. Before each level i, `root_clusters` contains the
    level i clusters formed from contractions in the previous level, along with any additional children of level
    i+1 del clusters that will get deleted. All root clusters have their parent pointer cleared already.
    The sequence `del_clusters` contains the level i+1 clusters that were parents of del clusters at the previous
    level, along with the parents of any nodes that started a new remove ancestor path due to a root cluster
    merging into a non-root cluster at the previous level. The sequence 'dir_updates' contains the list of the
    intial update edges that still exist in level i.
    (PHASE 1) First we update level i with 'dir_updates'. At level 0 this corresponds to updating the level 0
    forest. Since some clusters never get deleted, it is necessary to propagate this change to all levels.
    This must happen before the reclustering at level i, but after the del clusters are processed in phase 3,
    so that the recluster has accurate tree information, but phase 3 retains the structure of the previous tree.
    (PHASE 2) Next we, will determine the new set of contractions over the level i root clusters using an MIS. We
    then allocate new level i+1 parent nodes and update the parent pointers of the root clusters.
    These new parents become the `next_root_clusters` along with those determined phase 3. Any contraction
    with a non-root cluster adds the level i+1 parent to `del_clusters`.
    (PHASE 3) First we take the level i+1 `del_clusters`, and populate `next_del_clusters` by mapping them to their
    level i+2 parents. We also determine which ones will actually be deleted. During this step we delete the
    del clusters. We also find any level i+1 non-del children of the level i+2 clusters which will be deleted.
    These will become `next_root_clusters` at the next level i+1.
    (PHASE 4) In this phase we insert any edge incident to a level i root cluster into level i+1, if it wasn't
    contracted away. This fills the adjacency lists of our new parent clusters and also updates any level i+1
    clusters that remain from the previous tree. This must happen after phase 3 so that we do not modify the
    structure of the old tree before processing the del clusters. */
    while (root_clusters.size() > 0 || del_clusters.size() > 0) {
        // =======
        // PHASE 1
        // =======

        // For insertion batches, insert the linked edges at level i, and map to the next level.
        if (update_type == INSERT) {
            parlay::parallel_for(0, dir_update_groups.size(), [&] (size_t i) {
                auto& [cluster, edges] = dir_update_groups[i];
                auto neighbors = parlay::map(edges, [&] (auto x) { return x.second; });
                cluster->insert_neighbors(neighbors);
            });

            dir_updates = parlay::map_maybe(dir_updates, [&] (std::pair<Cluster*, Cluster*> x) -> std::optional<std::pair<Cluster*, Cluster*>> {
                if (x.first->parent && x.second->parent && x.first->parent != x.second->parent)
                    return std::make_pair(x.first->parent, x.second->parent);
                return std::nullopt;
            });
            dir_update_groups = integer_group_by_key_inplace(dir_updates);
        }

        // =======
        // PHASE 2
        // =======

        // Recluster the root clusters.
        auto additional_del_clusters = recluster_root_clusters(root_clusters);
        del_clusters.append(additional_del_clusters);

        // This returns only the new clusters that were created during the reclustering at this level.
        auto next_root_clusters_1 = create_new_parents(root_clusters);

        // =======
        // PHASE 3
        // =======

        // Group del clusters by parent.
        auto parents = parlay::tabulate(del_clusters.size(), [&] (size_t i) {
            return std::make_pair(del_clusters[i]->parent, del_clusters[i]);
        });
        auto parent_groups = integer_group_by_key_inplace(parents);

        // Get the next level i+2 del clusters.
        auto next_del_clusters = parlay::map_maybe(parent_groups, [&] (auto group) -> std::optional<Cluster*> {
            if (!group.first) return std::nullopt;
            return group.first;
        });

        // For deletion batches, compute new degrees of level i+1 del clusters.
        // For deletion batches, map the level i+1 updates to level i+2 before processing level i+1.
        if (update_type == DELETE) {
            parlay::parallel_for(0, dir_update_groups.size(), [&] (size_t i) {
                auto& [cluster, edges] = dir_update_groups[i];
                cluster->degree = cluster->get_degree() - edges.size();
            });

            temp_dir_updates = parlay::map_maybe(dir_updates, [&] (std::pair<Cluster*, Cluster*> x) -> std::optional<std::pair<Cluster*, Cluster*>> {
                if (x.first->parent && x.second->parent && x.first->parent != x.second->parent)
                    return std::make_pair(x.first->parent, x.second->parent);
                return std::nullopt;
            });
        }

        // Get the next level i+1 root clusters, which are children of a level i+2 del cluster that will be deleted.
        auto next_root_clusters_2 = process_del_clusters(parent_groups);
        parlay::parallel_for(0, next_root_clusters_2.size(), [&] (size_t i) {
            next_root_clusters_2[i]->parent = nullptr;
        });

        // For deletion batches, do the deletions at level i+1, finish mapping the updates to level i+2.
        if (update_type == DELETE) {
            parlay::parallel_for(0, dir_update_groups.size(), [&] (size_t i) {
                auto& [cluster, edges] = dir_update_groups[i];
                auto neighbors = parlay::map(edges, [&] (auto x) { return x.second; });
                cluster->delete_neighbors(neighbors);
            });

            dir_updates = std::move(temp_dir_updates);
            dir_update_groups = integer_group_by_key_inplace(dir_updates);
        }

        // Determine pointers to level i+1 del clusters that will get deleted.
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

        // Delete pointers to level i+1 del clusters that will get deleted.
        auto del_edge_groups = integer_group_by_key_inplace(del_edges);
        parlay::parallel_for(0, del_edge_groups.size(), [&] (size_t i) {
            auto& [cluster, edges] = del_edge_groups[i];
            auto neighbors = parlay::map(edges, [&] (auto x) { return x.second; });
            cluster->delete_neighbors(neighbors);
        });

        // Delete level i+1 del clusters that should be deleted.
        parlay::parallel_for(0, del_clusters.size(), [&] (size_t i) {
            if (del_clusters[i]->partner) allocator::destroy(del_clusters[i]);
        });

        // =======
        // PHASE 4
        // =======

        // Clear the partner fields of level i root cluster's and any partnered non-root clusters.
        // For all level i root clusters, insert each incident edge into level i+1 if possible.
        // This code uses the fine-grained locking insert and should work well for only low degree cases.
        parlay::parallel_for(0, root_clusters.size(), [&] (size_t i) {
            auto cluster = root_clusters[i];
            if (!cluster->parent) return; // Only deg 0

            // Clear partner pointers
            Cluster* partner = AtomicLoad(&cluster->partner);
            if (partner == (Cluster*) NEW_PAR_MARK) {
                AtomicStore(&cluster->partner, (Cluster*) NULL_PTR);
            }
            else if (partner) {
                Cluster* partner_partner = AtomicLoad(&partner->partner);
                if (partner_partner != cluster) { // Non-root partner
                    AtomicStore(&partner->partner, (Cluster*) NULL_PTR);
                    AtomicStore(&cluster->partner, (Cluster*) NULL_PTR);
                } else if (cluster < partner) { // Tie-break
                    AtomicStore(&partner->partner, (Cluster*) NULL_PTR);
                    AtomicStore(&cluster->partner, (Cluster*) NULL_PTR);
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

        // ==================
        // PREPARE NEXT LEVEL
        // ==================
        root_clusters = parlay::append(next_root_clusters_1, next_root_clusters_2);
        del_clusters = std::move(next_del_clusters);
    }
}

// ================================================================================================
// ================================        HELPER FUNCTIONS        ================================
// ================================================================================================

template <typename aug_t>
parlay::sequence<ParallelUFOCluster<aug_t>*> ParallelUFOTree<aug_t>::process_initial_clusters(parlay::sequence<std::pair<Cluster*, EdgeSlice>>& parent_groups) {
    return parlay::flatten(parlay::tabulate(parent_groups.size(), [&] (size_t i) {
        auto& [parent, children] = parent_groups[i];
        if (parent == nullptr) return parlay::map(children, [&] (auto child) { return child.second; });
        parlay::sequence<Cluster*> local_root_clusters;

        Cluster* max = (*parlay::max_element(children, [&] (auto x, auto y) { return x.second->get_degree() < y.second->get_degree(); })).second;
        size_t max_degree = max->get_degree();

        if (max_degree == 1) {
            Cluster* center = max->get_neighbor();
            if (center->parent != max->parent) {
                local_root_clusters.push_back(max);
                parent->partner = (Cluster*) DEL_MARK;
            } else if (center->degree < children.size() + 3) {
                local_root_clusters = center->filter_neighbors([&] (auto neighbor) {
                    return neighbor->parent == parent;
                });
                local_root_clusters.push_back(center);
                parent->partner = (Cluster*) DEL_MARK;
            } else {
                local_root_clusters = parlay::map(children, [&] (auto x) { return x.second; });
            }
        }

        else {
            if (max->degree <= (children.size()-1) + 3) {
                local_root_clusters = max->filter_neighbors([&] (auto neighbor) {
                    return neighbor->parent == parent;
                });
                local_root_clusters.push_back(max);
                parent->partner = (Cluster*) DEL_MARK;
            } else {
                local_root_clusters = parlay::map_maybe(children, [&] (auto child) -> std::optional<Cluster*> {
                    if (child.second != max) return child.second;
                    return std::nullopt;
                });
            }
        }

        return local_root_clusters;
    }));
}

template <typename aug_t>
parlay::sequence<ParallelUFOCluster<aug_t>*> ParallelUFOTree<aug_t>::process_del_clusters(parlay::sequence<std::pair<Cluster*, EdgeSlice>>& parent_groups) {
    return parlay::flatten(parlay::tabulate(parent_groups.size(), [&] (size_t i) {
        auto& [parent, children] = parent_groups[i];
        parlay::sequence<Cluster*> local_root_clusters;
        if (parent == nullptr) return parlay::map_maybe(children, [&] (auto child) -> std::optional<Cluster*> {
            if (!child.second->partner) return child.second;
            return std::nullopt;
        });

        Cluster* max = (*parlay::max_element(children, [&] (auto x, auto y) { return x.second->get_degree() < y.second->get_degree(); })).second;
        size_t max_degree = max->get_degree();

        if (max_degree == 1) {
            Cluster* center = max->get_neighbor();
            if (center->parent != max->parent) {
                if (!max->partner) local_root_clusters.push_back(max);
                parent->partner = (Cluster*) DEL_MARK;
            } else if (center->degree < children.size() + 3) {
                local_root_clusters = center->filter_neighbors([&] (auto neighbor) {
                    return !neighbor->partner && neighbor->parent == parent;
                });
                if (!center->partner) local_root_clusters.push_back(center);
                parent->partner = (Cluster*) DEL_MARK;
            } else {
                local_root_clusters = parlay::map_maybe(children, [&] (auto x) -> std::optional<Cluster*> {
                    if (!x.second->partner) return x.second;
                    return std::nullopt;
                });
            }
        }

        else {
            if (max->degree <= (children.size()-1) + 3) {
                local_root_clusters = max->filter_neighbors([&] (auto neighbor) {
                    return !neighbor->partner && neighbor->parent == parent;
                });
                if (!max->partner) local_root_clusters.push_back(max);
                parent->partner = (Cluster*) DEL_MARK;
            } else {
                local_root_clusters = parlay::map_maybe(children, [&] (auto x) -> std::optional<Cluster*> {
                    if (!x.second->partner && x.second != max) return x.second;
                    return std::nullopt;
                });
            }
        }

        return local_root_clusters;
    }));
}

template <typename aug_t>
parlay::sequence<ParallelUFOCluster<aug_t>*> ParallelUFOTree<aug_t>::recluster_root_clusters(parlay::sequence<Cluster*>& root_clusters) {
    // This function sets the partner fields for all root clusters. For a non-root
    // cluster, we mark its partner field with NON_ROOT_MARK. For root clusters that
    // don't combine with anything, we leave its partner field empty. For high degree
    // root clusters, we assign its partner field as NEW_PAR_MARK, and we add a parent
    // for it in this part. All other root clusters receive no parent at this point.
    // This returns the parent of any non-root clusters that were partnered with.
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
            neighbor->partner = (Cluster*) NON_ROOT_MARK;
            return neighbor->parent;
        }
    }
    else if (neighbor->get_degree() == 2) { // Combine deg 1 root cluster with deg 2 non-root clusters that don't contract
        if (neighbor->parent && !neighbor->contracts()) {
            if (CAS(&neighbor->partner, (Cluster*) NULL_PTR, (Cluster*) NON_ROOT_MARK)) {
                cluster->partner = neighbor;
                return neighbor->parent;
            }
        }
    }
    else { // Combine deg 1 root cluster with possible deg 3+ non-root clusters
        cluster->partner = neighbor;
        if (AtomicLoad(&neighbor->parent))
            if (CAS(&neighbor->partner, (Cluster*) NULL_PTR, (Cluster*) NON_ROOT_MARK))
                return neighbor->parent;
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
            if (curr != cluster && is_local_max(curr)) break;
            if (next->get_degree() == 2 && !next->parent && is_local_max(next)) break;
            if (!CAS(&curr->partner, (Cluster*) NULL_PTR, next)) break;
            if (next->get_degree() == 1) { // If next deg 1 they can definitely combine
                if (!next->parent) next->partner = curr;
                else next->partner = (Cluster*) NON_ROOT_MARK;
            } else { // deg 2
                Cluster* new_partner = next->parent ? (Cluster*) NON_ROOT_MARK : curr;
                if (!CAS(&next->partner, (Cluster*) NULL_PTR, new_partner)) { // If the CAS fails next was combined from the other side
                    if (AtomicLoad(&next->partner) != curr) // Other side combined with the opposite cluster (you got left hanging)
                        AtomicStore(&curr->partner, (Cluster*) NULL_PTR);
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
            else break;
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
    AtomicStore(&cluster->parent, parent);
    AtomicStore(&cluster->partner, (Cluster*) NEW_PAR_MARK);
    cluster->for_all_neighbors([&] (auto neighbor) {
        if (neighbor->get_degree() == 1 && neighbor->parent) {
            neighbor->parent->partner = (Cluster*) DEL_MARK;
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
        if (cluster->get_degree() >= 3 && cluster->partner == (Cluster*) NEW_PAR_MARK) return cluster->parent;
        Cluster* partner = cluster->partner;
        if (partner) {
            if (partner->partner != cluster) { // Non-root partner or high-degree partner with no partner field set
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
    if (neighbor1->get_degree() == 2 && !neighbor1->parent)
        if (hash1 > hash || (hash1 == hash && neighbor1 > c))
            return false;
    if (neighbor2->get_degree() == 2 && !neighbor2->parent)
        if (hash2 > hash || (hash2 == hash && neighbor2 > c))
            return false;
    return true;
}

template <typename aug_t>
bool ParallelUFOTree<aug_t>::connected(vertex_t u, vertex_t v) {
    return (leaves[u].get_root() == leaves[v].get_root());
}

}

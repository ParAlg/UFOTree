#pragma once
#include <parlay/alloc.h>
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include "bridge.h"
#include "sparse_table.h"
#include "types.h"
#include "util.h"


using namespace parlay;

// #define COLLECT_ROOT_CLUSTER_STATS
#ifdef COLLECT_ROOT_CLUSTER_STATS
std::map<int, int> root_clusters_histogram;
#endif
// #define COLLECT_HEIGHT_STATS
#ifdef COLLECT_HEIGHT_STATS
int max_height = 0;
#endif

long parallel_topology_remove_ancestor_time = 0;
long parallel_topology_remove_ancestor_time_2 = 0;
long parallel_topology_recluster_tree_time = 0;
long test_time = 0;

template <typename aug_t>
struct ParallelTopologyCluster {
  // Topology cluster data
  ParallelTopologyCluster<aug_t>* neighbors[3];
  ParallelTopologyCluster<aug_t>* parent;
  ParallelTopologyCluster<aug_t>* partner;

  aug_t edge_values[3];   // Only for path queries
  aug_t value;            // Stores subtree values or cluster path values
  // Some fields for asynchronous operations
  uint32_t priority;
  bool del;

  // Constructor
  ParallelTopologyCluster(){};
  ParallelTopologyCluster(uint32_t salt, aug_t value = 1)
      : neighbors(),
        edge_values(),
        value(value),
        parent(nullptr),
        del(false),
        partner(nullptr) {
    priority = hash32(salt);
  };
  // Helper functions
  int get_degree();
  bool contracts();
  bool contains_neighbor(ParallelTopologyCluster<aug_t>* c);
  void insert_neighbor(ParallelTopologyCluster<aug_t>* c, aug_t value);
  void remove_neighbor(ParallelTopologyCluster<aug_t>* c);
  ParallelTopologyCluster<aug_t>* get_neighbor(int index = 0);
  ParallelTopologyCluster<aug_t>* get_other_neighbor(
     ParallelTopologyCluster<aug_t>* first_neighbor);
  ParallelTopologyCluster<aug_t>* get_root();

  // Try to delete this cluster non-atomically.
  bool try_del();
  // Try to delete this cluster *atomically* (i.e., by CAS'ing the del
  // field).
  bool try_del_atomic();
};


template <typename aug_t>
class ParallelTopologyTree {
 public:

  using Cluster = ParallelTopologyCluster<aug_t>;
  using ClusterSubset = parlay::sequence<Cluster*>;

  // Topology tree interface
  ParallelTopologyTree(
     vertex_t n, vertex_t k, QueryType q = PATH,
     std::function<aug_t(aug_t, aug_t)> f = [](aug_t x, aug_t y) -> aug_t {
       return x + y;
     },
     aug_t id = 0, aug_t dval = 0);
  ~ParallelTopologyTree();
  void batch_link(sequence<Edge>& links);
  void batch_cut(sequence<Edge>& cuts);
  bool connected(vertex_t u, vertex_t v);
  aug_t subtree_query(vertex_t v, vertex_t p = MAX_VERTEX_T);
  aug_t path_query(vertex_t u, vertex_t v);
  // Testing helpers
  bool is_valid();
  int get_height(vertex_t v);
  void print_tree();

 private:
  // Class data and parameters
  vertex_t n;
  ParallelTopologyCluster<aug_t>* leaves;
  QueryType query_type;
  std::function<aug_t(aug_t, aug_t)> f;
  aug_t identity;
  aug_t default_value;
  parlay::sequence<ParallelTopologyCluster<aug_t>*> root_clusters;
  parlay::sequence<ParallelTopologyCluster<aug_t>*> del_clusters;
  // Helper functions
  void recluster_tree();
  void recompute_parent_value(ParallelTopologyCluster<aug_t>* c1,
                              ParallelTopologyCluster<aug_t>* c2);

  // Initializes the tree at the leaves.
  void initialize_from_leaves(sequence<Edge>& edges);
};

template <typename aug_t>
ParallelTopologyTree<aug_t>::ParallelTopologyTree(
   vertex_t n, vertex_t k, QueryType q, std::function<aug_t(aug_t, aug_t)> f,
   aug_t id, aug_t d)
    : n(n), query_type(q), f(f), identity(id), default_value(d) {
  leaves = new ParallelTopologyCluster<aug_t>[n];
  for (int i = 0; i < n; ++i)
    leaves[i] = ParallelTopologyCluster<aug_t>((uint32_t)i);
}

template <typename aug_t>
ParallelTopologyTree<aug_t>::~ParallelTopologyTree() {
  // Clear all memory
  std::unordered_set<ParallelTopologyCluster<aug_t>*> clusters;
  for (int i = 0; i < n; ++i) {
    auto curr = leaves[i].parent;
    while (curr) {
      clusters.insert(curr);
      curr = curr->parent;
    }
  }
  for (auto cluster : clusters)
    type_allocator<ParallelTopologyCluster<aug_t>>::destroy(cluster);
  delete[] leaves;
#ifdef COLLECT_ROOT_CLUSTER_STATS
  std::cout << "Number of root clusters: Frequency" << std::endl;
  for (auto entry : root_clusters_histogram)
    std::cout << entry.first << "\t" << entry.second << std::endl;
#endif
#ifdef COLLECT_HEIGHT_STATS
  std::cout << "Maximum height of the tree: " << max_height << std::endl;
#endif
  PRINT_TIMER("My Time", parallel_topology_remove_ancestor_time_2);
  PRINT_TIMER("REMOVE ANCESTORS TIME", parallel_topology_remove_ancestor_time);
  PRINT_TIMER("RECLUSTER TREE TIME", parallel_topology_recluster_tree_time);
  PRINT_TIMER("TEST TIME", test_time);
  return;
}

template <typename aug_t>
void ParallelTopologyTree<aug_t>::initialize_from_leaves(sequence<Edge>& edges) {
  START_TIMER(ra_timer);
  // Serial elision.
  if (edges.size() < 5000) {
    root_clusters.clear();
    del_clusters.clear();
    for (size_t i=0; i <2*edges.size(); ++i) {
      auto cluster = (i % 2 == 0) ? &leaves[edges[i / 2].src] : &leaves[edges[i / 2].dst];
      for (size_t j=0; j < 3; ++j) {
        auto neighbor = cluster->neighbors[j];
        if (neighbor && neighbor->parent == cluster->parent && !neighbor->del) {
          neighbor->del = true;
          root_clusters.push_back(neighbor);
          if (neighbor->parent && !neighbor->parent->del) {
            neighbor->parent->del = true;
            del_clusters.push_back(neighbor->parent);
          }
        }
      }
      if (!cluster->del) {
        cluster->del = true;
        root_clusters.push_back(cluster);
        if (cluster->parent && !cluster->parent->del) {
          cluster->parent->del = true;
          del_clusters.push_back(cluster->parent);
        }
      }
    }
  } else {
    auto all_clusters = parlay::delayed::tabulate(2*edges.size(), [&] (size_t i) {
      auto cluster = (i % 2 == 0) ? &leaves[edges[i / 2].src] : &leaves[edges[i / 2].dst];
      auto ds = parlay::delayed::tabulate(4, [cluster] (size_t j) {
        if (j < 3) {
          auto neighbor = cluster->neighbors[j];
          if (neighbor && neighbor->parent == cluster->parent) {
            return neighbor;
          }
          return (Cluster*)nullptr;
        } else {
          return cluster;
        }
      });
      return parlay::delayed::map_maybe(ds, [&] (auto cluster) -> std::optional<Cluster*>{
        if (cluster && cluster->try_del_atomic()) {
          return cluster;
        }
        return std::nullopt;
      });
    });
    root_clusters = parlay::delayed::to_sequence(parlay::delayed::flatten(all_clusters));

    del_clusters = map_maybe(
       root_clusters,
       [&](auto cluster) -> std::optional<ParallelTopologyCluster<aug_t>*> {
         if (cluster->parent && cluster->parent->try_del_atomic())
           return cluster->parent;
         return std::nullopt;
       });
  }
  STOP_TIMER(ra_timer, parallel_topology_remove_ancestor_time);
}


template <typename aug_t>
void ParallelTopologyTree<aug_t>::batch_link(sequence<Edge>& links) {
  // Generate root clusters.
  initialize_from_leaves(links);
  parallel_for(0, links.size(), [&](size_t i) {
    leaves[links[i].src].insert_neighbor(&leaves[links[i].dst], default_value);
    leaves[links[i].dst].insert_neighbor(&leaves[links[i].src], default_value);
  });
  recluster_tree();
}

template <typename aug_t>
void ParallelTopologyTree<aug_t>::batch_cut(sequence<Edge>& cuts) {
  initialize_from_leaves(cuts);
  parallel_for(0, cuts.size(), [&](size_t i) {
    leaves[cuts[i].src].remove_neighbor(&leaves[cuts[i].dst]);
    leaves[cuts[i].dst].remove_neighbor(&leaves[cuts[i].src]);
  });
  recluster_tree();
}

template <typename aug_t>
void ParallelTopologyTree<aug_t>::recluster_tree() {
  while (root_clusters.size() > 0 || del_clusters.size() > 0) {
// Update root cluster stats if we are collecting them
#ifdef COLLECT_ROOT_CLUSTER_STATS
    if (root_clusters_histogram.find(root_clusters.size()) ==
        root_clusters_histogram.end())
      root_clusters_histogram[root_clusters.size()] = 1;
    else
      root_clusters_histogram[root_clusters.size()] += 1;
#endif
    START_TIMER(ra_timer1);

    // Get the next set of clusters to delete
    parlay::sequence<Cluster*> new_del_clusters;
    parlay::sequence<Cluster*> next_additional_root_clusters;
    if (del_clusters.size() < 5000) {
      for (size_t i=0; i < del_clusters.size(); ++i) {
        auto cluster = del_clusters[i];
        if (cluster->parent && cluster->parent->try_del()) {
          new_del_clusters.push_back(cluster->parent);
        }
        if (cluster->parent) {
          for (auto neighbor : cluster->neighbors) {
            if (neighbor && !neighbor->del && neighbor->parent == cluster->parent && neighbor->try_del()) {
              next_additional_root_clusters.push_back(neighbor);
            }
          }
        }
      }
    } else {
      new_del_clusters = map_maybe(
         del_clusters,
         [&](auto cluster) -> std::optional<ParallelTopologyCluster<aug_t>*> {
           if (cluster->parent && cluster->parent->try_del_atomic())
             return cluster->parent;
           return std::nullopt;
         });
      // Determine the next level additional root clusters neighboring clusters we
      // are currently deleting
      next_additional_root_clusters = map_maybe(
         del_clusters,
         [&](auto cluster) -> std::optional<ParallelTopologyCluster<aug_t>*> {
           if (cluster->parent)
             for (auto neighbor : cluster->neighbors)
               if (neighbor && !neighbor->del &&
                   neighbor->parent == cluster->parent && neighbor->try_del_atomic())
                 return neighbor;
           return std::nullopt;
         });
      // Delete the current set of clusters to delete
    }

    // Remove links from nodes that are still live to nodes that will
    // be deleted.
    parallel_for(0, del_clusters.size(), [&](size_t i) {
      auto cluster = del_clusters[i];
      for (auto neighbor : cluster->neighbors)
        if (neighbor) {
          neighbor->remove_neighbor(cluster);
        }
    });
    // Destroy all nodes to be deleted.
    parallel_for(0, del_clusters.size(), [&](size_t i) {
      auto cluster = del_clusters[i];
      type_allocator<ParallelTopologyCluster<aug_t>>::destroy(cluster);
    });
    STOP_TIMER(ra_timer1, parallel_topology_remove_ancestor_time);

    START_TIMER(recluster_timer);
    // Set all root clusters parent to null
    parallel_for(0, root_clusters.size(), [&](size_t i) {
      root_clusters[i]->del = false;
      root_clusters[i]->parent = nullptr;
    });
    // Perform clustering of the root clusters
    parallel_for(0, root_clusters.size(), [&](size_t i) {
      auto cluster = root_clusters[i];
      // Combine deg 3 root clusters with deg 1 root  or non-root clusters
      if (cluster->get_degree() == 3) {
        for (auto neighbor : cluster->neighbors) {
          if (neighbor && neighbor->get_degree() == 1) {
            cluster->partner = neighbor;
            neighbor->partner = cluster;
            break;
          }
        }
      } else if (cluster->get_degree() == 2 && !cluster->parent) {
        // Only local maxima in priority with respect to deg 2 clusters will act
        bool local_max = true;
        for (auto neighbor : cluster->neighbors)
          if (neighbor && !neighbor->parent && neighbor->get_degree() == 2 &&
              neighbor->priority >= cluster->priority)
            local_max = false;
        if (!local_max)
          return;
        for (auto direction : {0, 1}) {
          // Travel left/right and pair clusters until a deg 3, deg 1, non-root,
          // or combined cluster found
          auto curr = cluster;
          auto next = cluster->get_neighbor(direction);
          if (curr->partner) {
            curr = cluster->get_neighbor(direction);
            next = curr->get_other_neighbor(cluster);
          }
          while (curr && !curr->parent && curr->get_degree() == 2 && next &&
                 next->get_degree() < 3 && !next->contracts()) {
            if (curr->partner ||
                !gbbs::CAS(&curr->partner, (ParallelTopologyCluster<aug_t>*)nullptr,
                     next))
              break;
            if (next->get_degree() == 1) {   // If next deg 1 they can combine
              next->partner = curr;
              break;
            }
            if (next->partner ||
                !gbbs::CAS(&next->partner, (ParallelTopologyCluster<aug_t>*)nullptr,
                     curr)) {   // If the CAS fails next was combined from the
                                // other side
              if (next->partner != curr)
                curr->partner = nullptr;
              break;
            }
            if (next->parent)
              break;   // Stop traversing at a non-root cluster
            // Get the next two clusters in the chain
            curr = next->get_other_neighbor(curr);
            if (curr)
              next = curr->get_other_neighbor(next);
          }
        }
      } else if (cluster->get_degree() == 1) {
        for (auto neighbor : cluster->neighbors) {
          // Combine deg 1 root clusters with deg 1 root or non-root clusters
          if (neighbor && neighbor->get_degree() == 1) {
            cluster->partner = neighbor;
            neighbor->partner = cluster;
            break;
          }
          // Combine deg 1 root cluster with deg 2 or 3 non-root clusters that
          // don't contract
          if (neighbor && neighbor->parent &&
              (neighbor->get_degree() == 2 || neighbor->get_degree() == 3)) {
            if (neighbor->contracts())
              continue;
            if (neighbor->partner ||
                !gbbs::CAS(&neighbor->partner,
                     (ParallelTopologyCluster<aug_t>*)nullptr, cluster))
              continue;
            cluster->partner = neighbor;
            neighbor->partner = cluster;
            break;
          }
        }
      }
    });
    // Create new clusters for each new combination
    parallel_for(0, root_clusters.size(), [&](size_t i) {
      auto cluster = root_clusters[i];
      // Make new parent clusters for all new combinations
      auto partner = cluster->partner;
      if (partner) {
        if (partner->parent &&
            partner->parent !=
               cluster->parent) {   // partner is non-root cluster
          cluster->parent = partner->parent;
          partner->partner = nullptr;
          recompute_parent_value(cluster, partner);
        } else if (cluster > partner) {   // higher address cluster will do the
                                          // combination
          auto parent = type_allocator<ParallelTopologyCluster<aug_t>>::create(
             cluster->priority + partner->priority, default_value);
          partner->parent = parent;
          cluster->parent = parent;
          recompute_parent_value(cluster, partner);
        }
      } else if (!cluster->parent &&
                 cluster->get_degree() >
                    0) {   // clusters that don't combine get a new parent
        auto parent = type_allocator<ParallelTopologyCluster<aug_t>>::create(
           cluster->priority, default_value);
        cluster->parent = parent;
        cluster->partner = cluster;
        parent->value = cluster->value;
      }
    });
    // Fill in the neighbor lists of the new clusters
    parallel_for(0, root_clusters.size(), [&](size_t i) {
      auto c1 = root_clusters[i];
      auto c2 = c1->partner;
      auto parent = c1->parent;
      if (!c2)
        return;
      bool new_parent = (c2->partner == c1);
      if (new_parent && c1 < c2)
        return;
      for (int i = 0; i < 3; ++i) {
        if (c1->neighbors[i] &&
            c1->neighbors[i] != c2) {   // Don't add c2's parent (itself)
          parent->insert_neighbor(c1->neighbors[i]->parent, c1->edge_values[i]);
          c1->neighbors[i]->parent->insert_neighbor(parent, c1->edge_values[i]);
        }
      }
      for (int i = 0; i < 3; ++i) {
        if (c2->neighbors[i] &&
            c2->neighbors[i] != c1) {   // Don't add c1's parent (itself)
          parent->insert_neighbor(c2->neighbors[i]->parent, c2->edge_values[i]);
          c2->neighbors[i]->parent->insert_neighbor(parent, c2->edge_values[i]);
        }
      }
      c2->partner = nullptr;
      c1->partner = nullptr;
    });
    STOP_TIMER(recluster_timer, parallel_topology_recluster_tree_time);

    // Finalize the new root clusters and new list of clusters to be deleted
    parlay::sequence<Cluster*> new_root_clusters;
      START_TIMER(ra_timer2);
    if (root_clusters.size() < 5000) {
      for (size_t i=0; i < root_clusters.size(); ++i) {
        auto cluster = root_clusters[i];
        if (cluster->parent && cluster->parent->try_del()) {
          auto new_cluster = cluster->parent;

          // Try to add the parent to the next set of del clusters.
          if (new_cluster->parent && new_cluster->parent->try_del()) {
            new_del_clusters.push_back(new_cluster->parent);
          }

          next_additional_root_clusters.push_back(new_cluster);
          for (auto neighbor : new_cluster->neighbors) {
            if (neighbor && neighbor->parent == new_cluster->parent && neighbor->try_del()) {
              next_additional_root_clusters.push_back(neighbor);
            }
          }
        }
      }
      root_clusters = std::move(next_additional_root_clusters);
      del_clusters = std::move(new_del_clusters);
    } else {
      auto new_root_clusters = map_maybe(
         root_clusters,
         [&](auto cluster) -> std::optional<ParallelTopologyCluster<aug_t>*> {
           if (cluster->parent && cluster->parent->try_del_atomic())
             return cluster->parent;
           return std::nullopt;
         });
      auto additional_root_clusters = map_maybe(
         new_root_clusters,
         [&](auto cluster) -> std::optional<ParallelTopologyCluster<aug_t>*> {
           if (cluster->parent)
             for (auto neighbor : cluster->neighbors)
               if (neighbor && neighbor->parent == cluster->parent && neighbor->try_del_atomic()) {
                 return neighbor;
               }
           return std::nullopt;
         });
      root_clusters = append(append(new_root_clusters, additional_root_clusters),
                next_additional_root_clusters);
      auto additional_del_clusters = map_maybe(
         new_root_clusters,
         [&](auto cluster) -> std::optional<ParallelTopologyCluster<aug_t>*> {
           if (cluster->parent && cluster->parent->try_del_atomic()) {
             return cluster->parent;
           }
           return std::nullopt;
         });
      del_clusters = append(new_del_clusters, additional_del_clusters);
    }
    STOP_TIMER(ra_timer2, parallel_topology_remove_ancestor_time_2);
  }
}

template <typename aug_t>
void ParallelTopologyTree<aug_t>::recompute_parent_value(
   ParallelTopologyCluster<aug_t>* c1, ParallelTopologyCluster<aug_t>* c2) {
  assert(c1->parent == c2->parent);
  auto parent = c1->parent;
  if (query_type == SUBTREE) {
    parent->value = f(c1->value, c2->value);
  } else if (query_type == PATH && c1->get_degree() == 2 &&
             c2->get_degree() == 2) {
    aug_t edge_val;
    for (int i = 0; i < 3; ++i)
      if (c1->neighbors[i] == c2)
        edge_val = c1->edge_values[i];
    parent->value = f(f(c1->value, c2->value), edge_val);
  }
}

template <typename aug_t>
int ParallelTopologyCluster<aug_t>::get_degree() {
  int deg = 0;
  for (auto neighbor : this->neighbors)
    if (neighbor)
      deg++;
  return deg;
}

// Helper function which returns whether this cluster combines with another
// cluster.
template <typename aug_t>
bool ParallelTopologyCluster<aug_t>::contracts() {
  if (!parent)
    return false;
  bool contracts = false;
  for (auto neighbor : this->neighbors)
    if (neighbor && neighbor->parent == this->parent)
      contracts = true;
  return contracts;
}

template <typename aug_t>
bool ParallelTopologyCluster<aug_t>::contains_neighbor(
   ParallelTopologyCluster<aug_t>* c) {
  for (int i = 0; i < 3; ++i)
    if (this->neighbors[i] == c)
      return true;
  return false;
}

template <typename aug_t>
void ParallelTopologyCluster<aug_t>::insert_neighbor(
   ParallelTopologyCluster<aug_t>* c, aug_t value) {
  if (this->contains_neighbor(c))
    return;
  for (int i = 0; i < 3; ++i) {
    if (gbbs::CAS(&this->neighbors[i], (ParallelTopologyCluster<aug_t>*)nullptr, c)) {
      this->edge_values[i] = value;
      return;
    } else if (this->neighbors[i] == c)
      return;
  }
  std::cerr << "No space to insert neighbor." << std::endl;
  std::abort();
}

template <typename aug_t>
void ParallelTopologyCluster<aug_t>::remove_neighbor(
   ParallelTopologyCluster<aug_t>* c) {
  for (int i = 0; i < 3; ++i) {
    if (this->neighbors[i] == c) {
      this->neighbors[i] = nullptr;
      return;
    }
  }
  std::cerr << "Neighbor to delete not found." << std::endl;
  std::abort();
}

template <typename aug_t>
ParallelTopologyCluster<aug_t>* ParallelTopologyCluster<aug_t>::get_neighbor(
   int index) {
  assert(get_degree() > index);
  int neighbors_seen = 0;
  for (auto neighbor : neighbors) {
    if (neighbor) {
      if (neighbors_seen == index)
        return neighbor;
      neighbors_seen++;
    }
  }
  return nullptr;
}

template <typename aug_t>
ParallelTopologyCluster<aug_t>*
ParallelTopologyCluster<aug_t>::get_other_neighbor(
   ParallelTopologyCluster<aug_t>* first_neighbor) {
  assert(contains_neighbor(first_neighbor));
  for (auto neighbor : neighbors)
    if (neighbor && neighbor != first_neighbor)
      return neighbor;
  return nullptr;
}

template <typename aug_t>
ParallelTopologyCluster<aug_t>* ParallelTopologyCluster<aug_t>::get_root() {
  ParallelTopologyCluster<aug_t>* curr = this;
  while (curr->parent)
    curr = curr->parent;
  return curr;
}

template <typename aug_t>
bool ParallelTopologyCluster<aug_t>::try_del() {
  if (!del) {
    del = true;
    return true;
  }
  return false;
}

template <typename aug_t>
bool ParallelTopologyCluster<aug_t>::try_del_atomic() {
  if (!del) {
    bool ret = gbbs::CAS(&del, false, true);
    return ret;
  }
  return false;
}

/* Return true if and only if there is a path from vertex u to
vertex v in the tree. */
template <typename aug_t>
bool ParallelTopologyTree<aug_t>::connected(vertex_t u, vertex_t v) {
  return leaves[u].get_root() == leaves[v].get_root();
}

/* Returns the value of the associative function f applied over
the augmented values for all the vertices in the subtree rooted
at v with respect to its parent p. If p = -1 (MAX_VERTEX_T) then
return the sum over the entire tree containing v. */
template <typename aug_t>
aug_t ParallelTopologyTree<aug_t>::subtree_query(vertex_t v, vertex_t p) {
  assert(v >= 0 && v < n && p >= 0 && (p < n || p == MAX_VERTEX_T));
  if (p == MAX_VERTEX_T)
    return leaves[v].get_root()->value;
  assert(leaves[v].contains_neighbor(&leaves[p]));
  // Get the total up until the LCA of v and p
  aug_t total = f(identity, leaves[v].value);
  auto curr_v = &leaves[v];
  auto curr_p = &leaves[p];
  while (curr_v->parent != curr_p->parent) {
    for (auto neighbor : curr_v->neighbors)
      if (neighbor && neighbor->parent == curr_v->parent)
        total =
           f(total, neighbor->value);   // Only count vertices on the side of v
    curr_v = curr_v->parent;
    curr_p = curr_p->parent;
  }
  // Add the total after the LCA of v and p
  if (curr_v->get_degree() == 1)
    return total;
  if (curr_v->get_degree() == 3) {
    curr_v = curr_v->parent;
    while (curr_v->parent) {
      for (auto neighbor : curr_v->neighbors)
        if (neighbor && neighbor->parent == curr_v->parent)
          total =
             f(total, neighbor->value);   // Count all remaining root clusters
      curr_v = curr_v->parent;
    }
    return total;
  }
  // If the cluster of v was deg 2 when it combined, only count the clusters on
  // the side of v
  ParallelTopologyCluster<aug_t>* curr_u;
  for (auto neighbor : curr_v->neighbors)
    if (neighbor && neighbor != curr_p)
      curr_u = neighbor;   // Find the neighbor of curr_v that is not curr_p
  while (curr_v->parent) {
    if (curr_u->parent == curr_v->parent) {
      total = f(
         total,
         curr_u->value);   // Count the remaining root clusters on the side of v
      if (curr_u->get_degree() == 1)
        return total;
      if (curr_u->get_degree() == 3) {
        curr_v = curr_v->parent;
        while (curr_v->parent) {
          for (auto neighbor : curr_v->neighbors)
            if (neighbor && neighbor->parent == curr_v->parent)
              total = f(total,
                        neighbor->value);   // Count all remaining root clusters
          curr_v = curr_v->parent;
        }
        return total;
      }
      for (auto neighbor : curr_u->neighbors)
        if (neighbor && neighbor != curr_v)
          curr_u =
             neighbor
                ->parent;   // Find the neighbor of curr_u that is not curr_v
      curr_v = curr_v->parent;
    } else {
      curr_u = curr_u->parent;
      curr_v = curr_v->parent;
    }
  }
  return total;
}

/* Returns the value of the associative function f applied over
the augmented values for all the edges on the unique path from
vertex u to vertex v. */
template <typename aug_t>
aug_t ParallelTopologyTree<aug_t>::path_query(vertex_t u, vertex_t v) {
  assert(u >= 0 && u < n && v >= 0 && v < n);
  assert(u != v && connected(u, v));
  // Compute the path on both sides for both vertices until they combine
  aug_t path_u1, path_u2, path_v1, path_v2;
  path_u1 = path_u2 = path_v1 = path_v2 = identity;
  ParallelTopologyCluster<aug_t>*bdry_u1, *bdry_u2, *bdry_v1, *bdry_v2;
  bdry_u1 = bdry_u2 = bdry_v1 = bdry_v2 = nullptr;
  if (leaves[u].get_degree() == 2) {
    bdry_u1 =
       leaves[u].neighbors[0] ? leaves[u].neighbors[0] : leaves[u].neighbors[1];
    bdry_u2 =
       leaves[u].neighbors[2] ? leaves[u].neighbors[2] : leaves[u].neighbors[1];
  }
  if (leaves[v].get_degree() == 2) {
    bdry_v1 =
       leaves[v].neighbors[0] ? leaves[v].neighbors[0] : leaves[v].neighbors[1];
    bdry_v2 =
       leaves[v].neighbors[2] ? leaves[v].neighbors[2] : leaves[v].neighbors[1];
  }
  auto curr_u = &leaves[u];
  auto curr_v = &leaves[v];
  while (curr_u->parent != curr_v->parent) {
    for (int i = 0; i < 3; ++i) {
      auto neighbor = curr_u->neighbors[i];
      if (neighbor && neighbor->parent == curr_u->parent) {
        if (curr_u->get_degree() == 2) {
          if (curr_u->parent->get_degree() == 2) {
            // Binary to Binary
            if (neighbor == bdry_u1) {
              path_u1 = f(path_u1, f(curr_u->edge_values[i], neighbor->value));
              bdry_u2 = bdry_u2->parent;
              for (int i = 0; i < 3; ++i)
                if (curr_u->parent->neighbors[i] &&
                    curr_u->parent->neighbors[i] != bdry_u2)
                  bdry_u1 = curr_u->parent->neighbors[i];
            } else {
              path_u2 = f(path_u2, f(curr_u->edge_values[i], neighbor->value));
              bdry_u1 = bdry_u1->parent;
              for (int i = 0; i < 3; ++i)
                if (curr_u->parent->neighbors[i] &&
                    curr_u->parent->neighbors[i] != bdry_u1)
                  bdry_u2 = curr_u->parent->neighbors[i];
            }
          } else {
            // Binary to Unary
            path_u1 = (neighbor == bdry_u1) ? path_u2 : path_u1;
          }
        } else {
          if (curr_u->parent->get_degree() == 2) {
            // Unary to Binary
            path_u1 = path_u2 = f(path_u1, curr_u->edge_values[i]);
            bdry_u1 = curr_u->parent->neighbors[0]
                         ? curr_u->parent->neighbors[0]
                         : curr_u->parent->neighbors[1];
            bdry_u2 = curr_u->parent->neighbors[2]
                         ? curr_u->parent->neighbors[2]
                         : curr_u->parent->neighbors[1];
          } else {
            // Unary to Unary
            path_u1 = f(path_u1, f(curr_u->edge_values[i], neighbor->value));
          }
        }
        break;
      }
    }
    if (!curr_u->contracts()) {
      if (bdry_u1)
        bdry_u1 = bdry_u1->parent;
      if (bdry_u2)
        bdry_u2 = bdry_u2->parent;
    }
    curr_u = curr_u->parent;
    for (int i = 0; i < 3; ++i) {
      auto neighbor = curr_v->neighbors[i];
      if (neighbor && neighbor->parent == curr_v->parent) {
        if (curr_v->get_degree() == 2) {
          if (curr_v->parent->get_degree() == 2) {
            // Binary to Binary
            if (neighbor == bdry_v1) {
              path_v1 = f(path_v1, f(curr_v->edge_values[i], neighbor->value));
              bdry_v2 = bdry_v2->parent;
              for (int i = 0; i < 3; ++i)
                if (curr_v->parent->neighbors[i] &&
                    curr_v->parent->neighbors[i] != bdry_v2)
                  bdry_v1 = curr_v->parent->neighbors[i];
            } else {
              path_v2 = f(path_v2, f(curr_v->edge_values[i], neighbor->value));
              bdry_v1 = bdry_v1->parent;
              for (int i = 0; i < 3; ++i)
                if (curr_v->parent->neighbors[i] &&
                    curr_v->parent->neighbors[i] != bdry_v1)
                  bdry_v2 = curr_v->parent->neighbors[i];
            }
          } else {
            // Binary to Unary
            path_v1 = (neighbor == bdry_v1) ? path_v2 : path_v1;
          }
        } else {
          if (curr_v->parent->get_degree() == 2) {
            // Unary to Binary
            path_v1 = path_v2 = f(path_v1, curr_v->edge_values[i]);
            bdry_v1 = curr_v->parent->neighbors[0]
                         ? curr_v->parent->neighbors[0]
                         : curr_v->parent->neighbors[1];
            bdry_v2 = curr_v->parent->neighbors[2]
                         ? curr_v->parent->neighbors[2]
                         : curr_v->parent->neighbors[1];
          } else {
            // Unary to Unary
            path_v1 = f(path_v1, f(curr_v->edge_values[i], neighbor->value));
          }
        }
        break;
      }
    }
    if (!curr_v->contracts()) {
      if (bdry_v1)
        bdry_v1 = bdry_v1->parent;
      if (bdry_v2)
        bdry_v2 = bdry_v2->parent;
    }
    curr_v = curr_v->parent;
  }
  // Get the correct path sides when the two vertices meet at the LCA
  aug_t total = identity;
  if (curr_u->get_degree() == 2)
    total = f(total, (curr_v == bdry_u1) ? path_u1 : path_u2);
  else
    total = f(total, path_u1);
  if (curr_v->get_degree() == 2)
    total = f(total, (curr_u == bdry_v1) ? path_v1 : path_v2);
  else
    total = f(total, path_v1);
  // Add the value of the last edge
  for (int i = 0; i < 3; ++i)
    if (curr_u->neighbors[i] == curr_v)
      total = f(total, curr_u->edge_values[i]);
  return total;
}

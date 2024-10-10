#pragma once
#include "types.h"
#include "util.h"
#include "ternarizable_interface.h"

// #define COLLECT_ROOT_CLUSTER_STATS
#ifdef COLLECT_ROOT_CLUSTER_STATS
    static std::map<int, int> root_clusters_histogram;
#endif
// #define COLLECT_HEIGHT_STATS
#ifdef COLLECT_HEIGHT_STATS
    static int max_height = 0;
#endif

static long topology_remove_ancestor_time = 0;
static long topology_recluster_tree_time = 0;

template<typename aug_t>
struct TopologyCluster {
    // Topology cluster data
    vertex_t neighbors[3];
    aug_t edge_values[3];   // Only for path queries
    aug_t value;            // Stores subtree values or cluster path values
    // Constructor
    TopologyCluster(aug_t value) : neighbors{NONE,NONE,NONE}, edge_values(), value(value){};
    // Helper functions
    int get_degree();
    bool contains_neighbor(vertex_t v);
    void insert_neighbor(vertex_t v, aug_t value);
    void remove_neighbor(vertex_t v);
};

template<typename aug_t>
class TopologyTree : ITernarizable {
public:
    // Topology tree interface
    TopologyTree(
     vertex_t n, QueryType q = PATH,
     std::function<aug_t(aug_t, aug_t)> f = [](aug_t x, aug_t y) -> aug_t {
       return x + y;}, aug_t id = 0, aug_t dval = 0);
    ~TopologyTree();
    void link(vertex_t u, vertex_t v, aug_t value);
    void link(vertex_t u, vertex_t v) { link(u,v,default_value); };
    void cut(vertex_t u, vertex_t v);
    void batch_link(Edge* links, int len);
    void batch_cut(Edge* cuts, int len);
    bool connected(vertex_t u, vertex_t v);
    aug_t subtree_query(vertex_t v, vertex_t p = MAX_VERTEX_T);
    aug_t path_query(vertex_t u, vertex_t v);
    
    //Interface methods overriden.
    short get_degree(vertex_t v) override {return clusters[v][0].get_degree();}
    std::pair<vertex_t, int> retrieve_v_to_del(vertex_t v) override {
      return std::pair(clusters[v][0].neighbors[0], clusters[v][0].edge_values[0]); 
    }
    // Testing helpers
    bool is_valid();
    int get_height(vertex_t v);
    void print_tree();
    size_t space();
private:
    // Class data and parameters
    std::vector<std::vector<TopologyCluster<aug_t>>> clusters;
    std::vector<vertex_t> parents;
    QueryType query_type;
    std::function<aug_t(aug_t, aug_t)> f;
    aug_t identity;
    aug_t default_value;
    std::vector<std::vector<vertex_t>> root_clusters;
    std::vector<std::pair<std::pair<vertex_t, vertex_t>, bool>> contractions;
    // Helper functions
    void remove_ancestors(vertex_t v, int start_level = 0);
    void recluster_tree();
    vertex_t get_parent(vertex_t v, int level);
    vertex_t get_root(vertex_t v);
    bool contracts(vertex_t v, int level);
    void recompute_parent_value(vertex_t c1, vertex_t c2, int level);
};

template<typename aug_t>
TopologyTree<aug_t>::TopologyTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d)
    : query_type(q), f(f), identity(id), default_value(d) {
    clusters.resize(n);
    for (int v = 0; v < n; ++v) clusters[v].resize(1, d);
    parents.resize(n, NONE);
    root_clusters.resize(max_tree_height(n));
    contractions.reserve(12);
}

template<typename aug_t>
TopologyTree<aug_t>::~TopologyTree() {
    #ifdef COLLECT_ROOT_CLUSTER_STATS
    std::cout << "Number of root clusters: Frequency" << std::endl;
        for (auto entry : root_clusters_histogram)
            std::cout << entry.first << "\t" << entry.second << std::endl;
    #endif
    #ifdef COLLECT_HEIGHT_STATS
        std::cout << "Maximum height of the tree: " << max_height << std::endl;
    #endif
    PRINT_TIMER("REMOVE ANCESTORS TIME", topology_remove_ancestor_time);
    PRINT_TIMER("RECLUSTER TREE TIME", topology_recluster_tree_time);
    return;
}

template<typename aug_t>
size_t TopologyTree<aug_t>::space(){ 
    std::unordered_set<TopologyCluster<aug_t>*> visited;
    size_t memory = sizeof(TopologyTree<aug_t>);
    for (auto leaf : clusters) {
        memory += sizeof(leaf);
        memory += leaf.size() * sizeof(TopologyCluster<aug_t>);
    }
    return memory;
}
/* Link vertex u and vertex v in the tree. Optionally include an
augmented value for the new edge (u,v). If no augmented value is
provided, the default value is 1. */
template<typename aug_t>
void TopologyTree<aug_t>::link(vertex_t u, vertex_t v, aug_t value) {
    assert(u >= 0 && u < clusters.size() && v >= 0 && v < clusters.size());
    assert(u != v && !connected(u,v));
    START_TIMER(topology_remove_ancestor_timer);
    remove_ancestors(u);
    remove_ancestors(v);
    STOP_TIMER(topology_remove_ancestor_timer, topology_remove_ancestor_time);
    clusters[u][0].insert_neighbor(v, value);
    clusters[v][0].insert_neighbor(u, value);
    START_TIMER(topology_recluster_tree_timer);
    recluster_tree();
    STOP_TIMER(topology_recluster_tree_timer, topology_recluster_tree_time);
    // Collect tree height stats at the end of each update
    #ifdef COLLECT_HEIGHT_STATS
        max_height = std::max(max_height, get_height(u));
        max_height = std::max(max_height, get_height(v));
    #endif
}

/* Cut vertex u and vertex v in the tree. */
template<typename aug_t>
void TopologyTree<aug_t>::cut(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < clusters.size() && v >= 0 && v < clusters.size());
    assert(clusters[u][0].contains_neighbor(v));
    START_TIMER(topology_remove_ancestor_timer);
    remove_ancestors(u);
    remove_ancestors(v);
    STOP_TIMER(topology_remove_ancestor_timer, topology_remove_ancestor_time);
    clusters[u][0].remove_neighbor(v);
    clusters[v][0].remove_neighbor(u);
    START_TIMER(topology_recluster_tree_timer);
    recluster_tree();
    STOP_TIMER(topology_recluster_tree_timer, topology_recluster_tree_time);
    // Collect tree height stats at the end of each update
    #ifdef COLLECT_HEIGHT_STATS
        max_height = std::max(max_height, get_height(u));
        max_height = std::max(max_height, get_height(v));
    #endif
}

template<typename aug_t>
void TopologyTree<aug_t>::batch_link(Edge* links, int len) {
    START_TIMER(topology_remove_ancestor_timer);
    for (int i = 0; i < len; i++) {
        Edge e = links[i];
        vertex_t u = e.src;
        vertex_t v = e.dst;
        remove_ancestors(u);
        remove_ancestors(v);
        clusters[u][0].insert_neighbor(v, default_value);
        clusters[v][0].insert_neighbor(u, default_value);
    }
    STOP_TIMER(topology_remove_ancestor_timer, topology_remove_ancestor_time);
    START_TIMER(topology_recluster_tree_timer);
    recluster_tree();
    STOP_TIMER(topology_recluster_tree_timer, topology_recluster_tree_time);
}

template<typename aug_t>
void TopologyTree<aug_t>::batch_cut(Edge* cuts, int len) {
    START_TIMER(topology_remove_ancestor_timer);
    for (int i = 0; i < len; i++) {
        Edge e = cuts[i];
        vertex_t u = e.src;
        vertex_t v = e.dst;
        remove_ancestors(u);
        remove_ancestors(v);
        clusters[u][0].remove_neighbor(v);
        clusters[v][0].remove_neighbor(u);
    }
    STOP_TIMER(topology_remove_ancestor_timer, topology_remove_ancestor_time);
    START_TIMER(topology_recluster_tree_timer);
    recluster_tree();
    STOP_TIMER(topology_recluster_tree_timer, topology_recluster_tree_time);
}

template<typename aug_t>
void TopologyTree<aug_t>::remove_ancestors(vertex_t v, int start_level) {
    int level = start_level;
    vertex_t curr = v;
    vertex_t p = get_parent(curr, level);
    for (auto neighbor : clusters[curr][level].neighbors) {
        if (neighbor != NONE) {
            auto neighbor_p = get_parent(neighbor, level);
            if (neighbor_p == p) {
                if (curr == p)
                    parents[neighbor] = NONE; // Set sibling parent pointer to NONE
                else
                    parents[curr] = NONE;
                root_clusters[level].push_back(neighbor);
            }
        }
    }
    root_clusters[level].push_back(curr);
    curr = p;
    int del_level = level+1;
    while (curr != NONE) {
        level++;
        vertex_t prev = curr;
        curr = get_parent(prev, level);
        for (auto neighbor : clusters[prev][level].neighbors) {
            if (neighbor != NONE) {
                auto neighbor_p = get_parent(neighbor, level);
                if (neighbor_p == curr) {
                    if (prev == curr)
                        parents[neighbor] = NONE; // Set sibling parent pointer to NONE
                    else
                        parents[prev] = NONE;
                    root_clusters[level].push_back(neighbor);
                }
                clusters[neighbor][level].remove_neighbor(prev);   // Remove prev from adjacency
            }
        }
        auto position = std::find(root_clusters[level].begin(), root_clusters[level].end(), prev);
        if (position != root_clusters[level].end())
            root_clusters[level].erase(position);
        if (curr != prev) { // delete some clusters
            assert(del_level < clusters[prev].size());
            clusters[prev].resize(del_level, default_value);
            del_level = level+1;
        }
    }
}

template<typename aug_t>
void TopologyTree<aug_t>::recluster_tree() {
    for (int level = 0; level < root_clusters.size(); level++) {
    if (root_clusters[level].empty())
      continue;
        // Update root cluster stats if we are collecting them
        #ifdef COLLECT_ROOT_CLUSTER_STATS
    if (root_clusters_histogram.find(root_clusters[level].size()) ==
        root_clusters_histogram.end())
                root_clusters_histogram[root_clusters[level].size()] = 1;
            else
                root_clusters_histogram[root_clusters[level].size()] += 1;
        #endif
        for (auto cluster_id : root_clusters[level]) {
            auto cluster = &clusters[cluster_id][level];
            vertex_t p_id = get_parent(cluster_id, level);
            if (cluster->get_degree() == 3 && p_id == NONE) {
                // Combine deg 3 root clusters with deg 1 root or non-root clusters
                for (auto neighbor_id : cluster->neighbors) {
                    if (neighbor_id == NONE) continue;
                    auto neighbor = &clusters[neighbor_id][level];
                    if (neighbor->get_degree() == 1) {
                        auto parent_id = get_parent(neighbor_id, level);
                        bool new_parent = (parent_id == NONE);
                        if (new_parent) { // If neighbor is a root cluster
                            parent_id = cluster_id;
                            clusters[cluster_id].push_back(default_value);
                            root_clusters[level+1].push_back(parent_id);
                            parents[neighbor_id] = parent_id;
                        } else parents[cluster_id] = parent_id;
                        recompute_parent_value(cluster_id, neighbor_id, level);
                        contractions.push_back({{cluster_id, neighbor_id}, new_parent});
                        break;
                    }
                }
            } else if (cluster->get_degree() == 2 && p_id == NONE) {
                // Combine deg 2 root clusters with deg 1 or 2 root clusters
                for (auto neighbor_id : cluster->neighbors) {
                    if (neighbor_id == NONE) continue;
                    auto neighbor = &clusters[neighbor_id][level];
                    auto neighbor_p_id = get_parent(neighbor_id, level);
                    if (neighbor_p_id == NONE && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
                        auto parent_id = cluster_id;
                        clusters[cluster_id].push_back(default_value);
                        parents[neighbor_id] = parent_id;
                        root_clusters[level+1].push_back(parent_id);
                        recompute_parent_value(cluster_id, neighbor_id, level);
                        contractions.push_back({{cluster_id, neighbor_id}, true});
                        p_id = parent_id;
                        break;
                    }
                }
                // Combine deg 2 root clusters with deg 1 or 2 non-root clusters
                if (p_id == NONE)
                for (auto neighbor_id : cluster->neighbors) {
                    if (neighbor_id == NONE) continue;
                    auto neighbor = &clusters[neighbor_id][level];
                    auto neighbor_p_id = get_parent(neighbor_id, level);
                    if (neighbor_p_id != NONE && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
                        if (contracts(neighbor_id, level)) continue;
                        parents[cluster_id] = neighbor_p_id;
                        recompute_parent_value(cluster_id, neighbor_id, level);
                        contractions.push_back({{cluster_id, neighbor_id}, false});
                        break;
                    }
                }
            } else if (cluster->get_degree() == 1 && p_id == NONE) {
                // Combine deg 1 root clusters with deg 1, 2, or 3 root clusters
                for (auto neighbor_id : cluster->neighbors) {
                    if (neighbor_id == NONE) continue;
                    auto neighbor = &clusters[neighbor_id][level];
                    auto neighbor_p_id = get_parent(neighbor_id, level);
                    if (neighbor_p_id == NONE) {
                        auto parent_id = cluster_id;
                        clusters[cluster_id].push_back(default_value);
                        parents[neighbor_id] = parent_id;
                        root_clusters[level+1].push_back(parent_id);
                        recompute_parent_value(cluster_id, neighbor_id, level);
                        contractions.push_back({{cluster_id, neighbor_id}, true});
                        p_id = parent_id;
                        break;
                    }
                }
                // Combine deg 1 root clusters with deg 1, 2, or 3 non-root clusters
                if (p_id == NONE)
                for (auto neighbor_id : cluster->neighbors) {
                    if (neighbor_id == NONE) continue;
                    auto neighbor = &clusters[neighbor_id][level];
                    auto neighbor_p_id = get_parent(neighbor_id, level);
                    if (neighbor_p_id != NONE) {
                        if (contracts(neighbor_id, level)) continue;
                        parents[cluster_id] = neighbor_p_id;
                        recompute_parent_value(cluster_id, neighbor_id, level);
                        contractions.push_back({{cluster_id, neighbor_id}, false});
                        break;
                    }
                }
            }
            // Add remaining uncombined clusters to the next level
            p_id = get_parent(cluster_id, level);
            if (p_id == NONE && cluster->get_degree() > 0) {
                clusters[cluster_id].push_back(cluster->value);
                root_clusters[level+1].push_back(cluster_id);
                contractions.push_back({{cluster_id, cluster_id}, true});
            }
        }
        // Fill in the neighbor lists of the new clusters
        for (auto contraction : contractions) {
            auto c1_id = contraction.first.first;
            auto c1 = &clusters[c1_id][level];
            auto c2_id = contraction.first.second;
            auto c2 = &clusters[c2_id][level];
            auto parent_id = get_parent(c1_id, level);
            auto parent = &clusters[parent_id][level+1];
            bool new_parent = contraction.second;
            for (int i = 0; i < 3; ++i) parent->neighbors[i] = NONE;
            for (int i = 0; i < 3; ++i) {
                if (c1->neighbors[i] != NONE && c1->neighbors[i] != c2_id) {   // Don't add c2's parent (itself)
                    auto neighbor_p_id = get_parent(c1->neighbors[i], level);
                    auto neighbor_p = &clusters[neighbor_p_id][level+1];
                    parent->insert_neighbor(neighbor_p_id, c1->edge_values[i]);
                    neighbor_p->insert_neighbor(parent_id, c1->edge_values[i]);
                }
            }
            for (int i = 0; i < 3; ++i) {
                if (c2->neighbors[i] != NONE && c2->neighbors[i] != c1_id) {   // Don't add c1's parent (itself)
                    auto neighbor_p_id = get_parent(c2->neighbors[i], level);
                    auto neighbor_p = &clusters[neighbor_p_id][level+1];
                    parent->insert_neighbor(neighbor_p_id, c2->edge_values[i]);
                    neighbor_p->insert_neighbor(parent_id, c2->edge_values[i]);
                }
            }
            if (!new_parent) remove_ancestors(parent_id, level+1);
        }
        // Clear the contents of this level
        root_clusters[level].clear();
        contractions.clear();
    }
}

template<typename aug_t>
vertex_t TopologyTree<aug_t>::get_parent(vertex_t v, int level) {
    if (level+1 < clusters[v].size()) return v;
    return parents[v];
}

template<typename aug_t>
vertex_t TopologyTree<aug_t>::get_root(vertex_t v) {
    vertex_t curr = v;
    while (parents[curr] != NONE) curr = parents[curr];
    return curr;
}

// Helper function which returns whether this cluster combines with another cluster
template<typename aug_t>
bool TopologyTree<aug_t>::contracts(vertex_t v, int level) {
    auto cluster_p = get_parent(v, level);
    bool contracts = false;
    for (auto neighbor : clusters[v][level].neighbors) {
        if (neighbor == NONE) continue;
        auto neighbor_p = get_parent(neighbor, level);
        if (neighbor_p == cluster_p) contracts = true;
    }
    return contracts;
}

template<typename aug_t>
void TopologyTree<aug_t>::recompute_parent_value(vertex_t c1, vertex_t c2, int level) {
    assert(get_parent(c1, level) == get_parent(c2, level));
    auto cluster1 = &clusters[c1][level];
    auto cluster2 = &clusters[c2][level];
    auto parent = &clusters[get_parent(c1, level)][level+1];
    if (query_type == SUBTREE) {
        parent->value = f(cluster1->value, cluster2->value);
    } else if (query_type == PATH && cluster1->get_degree() == 2 && cluster2->get_degree() == 2) {
        aug_t edge_val;
        for (int i = 0; i < 3; i++) if (cluster1->neighbors[i] == c2) edge_val = cluster1->edge_values[i];
        parent->value = f(f(cluster1->value, cluster2->value), edge_val);
    }
}

template<typename aug_t>
int TopologyCluster<aug_t>::get_degree() {
    int deg = 0;
    for (auto neighbor : this->neighbors)
        if (neighbor != NONE) deg++;
    return deg;
}

template<typename aug_t>
bool TopologyCluster<aug_t>::contains_neighbor(vertex_t v) {
    for (int i = 0; i < 3; ++i)
        if (this->neighbors[i] == v)
            return true;
    return false;
}

template<typename aug_t>
void TopologyCluster<aug_t>::insert_neighbor(vertex_t v, aug_t value) {
    if (this->contains_neighbor(v)) return;
    for (int i = 0; i < 3; ++i) {
        if (this->neighbors[i] == NONE) {
            this->neighbors[i] = v;
            this->edge_values[i] = value;
            return;
        }
    }
    throw std::invalid_argument("No space to insert neighbor"); 
}

template<typename aug_t>
void TopologyCluster<aug_t>::remove_neighbor(vertex_t v) {
    for (int i = 0; i < 3; ++i) {
        if (this->neighbors[i] == v) {
            this->neighbors[i] = NONE;
            return;
        }
    }
    throw std::invalid_argument("Neighbor to delete not found");
}

/* Return true if and only if there is a path from vertex u to
vertex v in the tree. */
template<typename aug_t>
bool TopologyTree<aug_t>::connected(vertex_t u, vertex_t v) {
    return get_root(u) == get_root(v);
}

/* Returns the value of the associative function f applied over
the augmented values for all the vertices in the subtree rooted
at v with respect to its parent p. If p = -1 (MAX_VERTEX_T) then
return the sum over the entire tree containing v. */
template<typename aug_t>
aug_t TopologyTree<aug_t>::subtree_query(vertex_t v, vertex_t p) {
    assert(v >= 0 && v < clusters.size() && p >= 0 && (p < clusters.size() || p == MAX_VERTEX_T));
    if (p == MAX_VERTEX_T) return clusters[get_root(v)].back().value;
    assert(clusters[v][0].contains_neighbor(p));
    // Get the total up until the LCA of v and p
    aug_t total = f(identity, clusters[v][0].value);
    auto curr_v = v;
    auto curr_p = p;
    int level = 0;
    auto parent_v = get_parent(curr_v, level);
    auto parent_p = get_parent(curr_p, level);
    while (parent_v != parent_p) {
        for (auto neighbor : clusters[curr_v][level].neighbors)
            if (neighbor != NONE && get_parent(neighbor, level) == parent_v)
                total = f(total, clusters[neighbor][level].value);   // Only count vertices on the side of v
        curr_v = parent_v;
        curr_p = parent_p;
        level++;
        parent_v = get_parent(curr_v, level);
        parent_p = get_parent(curr_p, level);
    }
    // Add the total after the LCA of v and p
    if (clusters[curr_v][level].get_degree() == 1) return total;
    if (clusters[curr_v][level].get_degree() == 3) {
        curr_v = parent_v;
        level++;
        parent_v = get_parent(parent_v, level);
        while (parent_v != NONE) {
            for (auto neighbor : clusters[curr_v][level].neighbors)
                if (neighbor != NONE && get_parent(neighbor, level) == parent_v)
                    total = f(total, clusters[neighbor][level].value);   // Count all remaining root clusters
            curr_v = parent_v;
            level++;
            parent_v = get_parent(curr_v, level);
        }
        return total;
    }
    // If the cluster of v was deg 2 when it combined, only count the clusters on the side of v
    vertex_t curr_u;
    for (auto neighbor : clusters[curr_v][level].neighbors)
        if (neighbor != NONE && neighbor != curr_p)
            curr_u = neighbor; // Find the neighbor of curr_v that is not curr_p
    vertex_t parent_u = get_parent(curr_u, level);
    while (parent_v != NONE) {
        if (parent_u == parent_v) {
            total = f(total, clusters[curr_u][level].value); // Count remaining root clusters on the side of v
            if (clusters[curr_u][level].get_degree() == 1) return total;
            if (clusters[curr_u][level].get_degree() == 3) {
                curr_v = parent_v;
                level++;
                parent_v = get_parent(parent_v, level);
                while (parent_v != NONE) {
                    for (auto neighbor : clusters[curr_v][level].neighbors)
                        if (neighbor != NONE && get_parent(neighbor, level) == parent_v)
                            total = f(total, clusters[neighbor][level].value);   // Count all remaining root clusters
                    curr_v = parent_v;
                    level++;
                    parent_v = get_parent(curr_v, level);
                }
                return total;
            }
            for (auto neighbor : clusters[curr_u][level].neighbors)
                if (neighbor != NONE && neighbor != curr_v)
                    curr_u = neighbor; // Find the neighbor of curr_u that is not curr_v
            curr_u = get_parent(curr_u, level);
            curr_v = parent_v;
            level++;
            parent_u = get_parent(curr_u, level);
            parent_v = get_parent(curr_v, level);
        } else {
            curr_u = parent_u;
            curr_v = parent_v;
            level++;
            parent_u = get_parent(curr_u, level);
            parent_v = get_parent(curr_v, level);
        }
    }
    return total;
}


/* Returns the value of the associative function f applied over
the augmented values for all the edges on the unique path from
vertex u to vertex v. */
template<typename aug_t>
aug_t TopologyTree<aug_t>::path_query(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < clusters.size() && v >= 0 && v < clusters.size());
    assert(u != v && connected(u,v));
    // Compute the path on both sides for both vertices until they combine
    aug_t path_u1, path_u2, path_v1, path_v2;
    path_u1 = path_u2 = path_v1 = path_v2 = identity;
    vertex_t bdry_u1, bdry_u2, bdry_v1, bdry_v2;
    bdry_u1 = bdry_u2 = bdry_v1 = bdry_v2 = NONE;
    if (clusters[u][0].get_degree() == 2) {
        bdry_u1 = clusters[u][0].neighbors[0] != NONE ? clusters[u][0].neighbors[0] : clusters[u][0].neighbors[1];
        bdry_u2 = clusters[u][0].neighbors[2] != NONE ? clusters[u][0].neighbors[2] : clusters[u][0].neighbors[1];
    }
    if (clusters[v][0].get_degree() == 2) {
        bdry_v1 = clusters[v][0].neighbors[0] != NONE ? clusters[v][0].neighbors[0] : clusters[v][0].neighbors[1];
        bdry_v2 = clusters[v][0].neighbors[2] != NONE ? clusters[v][0].neighbors[2] : clusters[v][0].neighbors[1];
    }
    auto curr_u = u;
    auto curr_v = v;
    int level = 0;
    auto parent_u = get_parent(curr_u, level);
    auto parent_v = get_parent(curr_v, level);
    while (parent_u != parent_v) { 
        // NOTE(ATHARVA): Make this all into one function.
        for (int i = 0; i < 3; i++) {
            auto neighbor = clusters[curr_u][level].neighbors[i];
            if (neighbor == NONE) continue;
            auto parent_n = get_parent(neighbor, level);
            if (parent_n == parent_u) {
                if (clusters[curr_u][level].get_degree() == 2) {
                    if (clusters[parent_u][level+1].get_degree() == 2) {
                        // Binary to Binary
                        if (neighbor == bdry_u1) {
                            path_u1 = f(path_u1, f(clusters[curr_u][level].edge_values[i], clusters[neighbor][level].value));
                            bdry_u2 = get_parent(bdry_u2, level);
                            for (auto parent_neighbor : clusters[parent_u][level+1].neighbors)
                                if (parent_neighbor != NONE && parent_neighbor != bdry_u2)
                                    bdry_u1 = parent_neighbor;
                        } else {
                            path_u2 = f(path_u2, f(clusters[curr_u][level].edge_values[i], clusters[neighbor][level].value));
                            bdry_u1 = get_parent(bdry_u1, level);
                            for (auto parent_neighbor : clusters[parent_u][level+1].neighbors)
                                if (parent_neighbor != NONE && parent_neighbor != bdry_u1)
                                    bdry_u2 = parent_neighbor;
                        }
                    } else {
                        // Binary to Unary
                        path_u1 = (neighbor == bdry_u1) ? path_u2 : path_u1;
                    }
                } else {
                    if (clusters[parent_u][level+1].get_degree() == 2) {
                        // Unary to Binary and degree 3 case
                        if(clusters[curr_u][level].get_degree() != 3) path_u1 = path_u2 = f(path_u1, clusters[curr_u][level].edge_values[i]);
                        bdry_u1 = clusters[parent_u][level+1].neighbors[0] ? clusters[parent_u][level+1].neighbors[0] : clusters[parent_u][level+1].neighbors[1];
                        bdry_u2 = clusters[parent_u][level+1].neighbors[2] ? clusters[parent_u][level+1].neighbors[2] : clusters[parent_u][level+1].neighbors[1];
                    } else {
                        // Unary to Unary
                        path_u1 = f(path_u1, f(clusters[curr_u][level].edge_values[i], clusters[neighbor][level].value));
                    }
                }
                break;
            }
        }
        for (int i = 0; i < 3; i++) {
            auto neighbor = clusters[curr_v][level].neighbors[i];
            if (neighbor == NONE) continue;
            auto parent_n = get_parent(neighbor, level);
            if (parent_n == parent_v) {
                if (clusters[curr_v][level].get_degree() == 2) {
                    if (clusters[parent_v][level+1].get_degree() == 2) {
                        // Binary to Binary
                        if (neighbor == bdry_v1) {
                            path_v1 = f(path_v1, f(clusters[curr_v][level].edge_values[i], clusters[neighbor][level].value));
                            bdry_v2 = get_parent(bdry_v2, level);
                            for (auto parent_neighbor : clusters[parent_v][level+1].neighbors)
                                if (parent_neighbor != NONE && parent_neighbor != bdry_v2)
                                    bdry_v1 = parent_neighbor;
                        } else {
                            path_v2 = f(path_v2, f(clusters[curr_v][level].edge_values[i], clusters[neighbor][level].value));
                            bdry_v1 = get_parent(bdry_v1, level);
                            for (auto parent_neighbor : clusters[parent_v][level+1].neighbors)
                                if (parent_neighbor != NONE && parent_neighbor != bdry_v1)
                                    bdry_v2 = parent_neighbor;
                        }
                    } else {
                        // Binary to Unary
                        path_v1 = (neighbor == bdry_v1) ? path_v2 : path_v1;
                    }
                } else {
                    if (clusters[parent_v][level+1].get_degree() == 2) {
                        // Unary to Binary and degree 3 case
                        if(clusters[curr_v][level].get_degree() != 3) path_v1 = path_v2 = f(path_v1, clusters[curr_v][level].edge_values[i]);
                        bdry_v1 = clusters[parent_v][level+1].neighbors[0] ? clusters[parent_v][level+1].neighbors[0] : clusters[parent_v][level+1].neighbors[1];
                        bdry_v2 = clusters[parent_v][level+1].neighbors[2] ? clusters[parent_v][level+1].neighbors[2] : clusters[parent_v][level+1].neighbors[1];
                    } else {
                        // Unary to Unary
                        path_v1 = f(path_v1, f(clusters[curr_v][level].edge_values[i], clusters[neighbor][level].value));
                    }
                }
                break;
            }
        }
        if (!contracts(curr_u, level)) {
            if (bdry_u1 != NONE) bdry_u1 = get_parent(bdry_u1, level);
            if (bdry_u2 != NONE) bdry_u2 = get_parent(bdry_u2, level);
        }
        if (!contracts(curr_v, level)) {
            if (bdry_v1 != NONE) bdry_v1 = get_parent(bdry_v1, level);
            if (bdry_v2 != NONE) bdry_v2 = get_parent(bdry_v2, level);
        }
        curr_u = parent_u;
        curr_v = parent_v;
        level++;
        parent_u = get_parent(curr_u, level);
        parent_v = get_parent(curr_v, level);
    }
    // Get the correct path sides when the two vertices meet at the LCA
    aug_t total = identity;
    if (clusters[curr_u][level].get_degree() == 2)
        total = f(total, (curr_v == bdry_u1) ? path_u1 : path_u2);
    else
        total = f(total, path_u1);
    if (clusters[curr_v][level].get_degree() == 2)
        total = f(total, (curr_u == bdry_v1) ? path_v1 : path_v2);
    else
        total = f(total, path_v1);
    // Add the value of the last edge
    for (int i = 0; i < 3; i++)
        if (clusters[curr_u][level].neighbors[i] == curr_v)
            total = f(total, clusters[curr_u][level].edge_values[i]);
    return total;
}

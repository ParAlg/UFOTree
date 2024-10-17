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

template<typename v_t, typename e_t>
struct TopologyCluster {
    // Topology cluster data
    TopologyCluster<v_t, e_t>* parent;
    TopologyCluster<v_t, e_t>* neighbors[3];
    e_t edge_values[3];   // Only for path queries
    v_t value;            // Stores subtree values or cluster path values
    // Constructor
    TopologyCluster() : parent(), neighbors(), edge_values(), value(){};
    // Helper functions
    int get_degree();
    bool contracts();
    bool contains_neighbor(TopologyCluster<v_t, e_t>* c);
    void insert_neighbor(TopologyCluster<v_t, e_t>* c, e_t value);
    void remove_neighbor(TopologyCluster<v_t, e_t>* c);
    TopologyCluster<v_t, e_t>* get_root();
};

template<typename v_t>
struct TopologyCluster<v_t, empty_t> {
    // Topology cluster data
    TopologyCluster<v_t, empty_t>* parent;
    TopologyCluster<v_t, empty_t>* neighbors[3];
    v_t value;            // Stores subtree values or cluster path values
    // Constructor
    TopologyCluster() : parent(), neighbors(), value(){};
    // Helper functions
    int get_degree();
    bool contracts();
    bool contains_neighbor(TopologyCluster<v_t, empty_t>* c);
    void insert_neighbor(TopologyCluster<v_t, empty_t>* c);
    void remove_neighbor(TopologyCluster<v_t, empty_t>* c);
    TopologyCluster<v_t, empty_t>* get_root();
};

template<>
struct TopologyCluster<empty_t, empty_t> {
    // Topology cluster data
    TopologyCluster<empty_t, empty_t>* parent;
    TopologyCluster<empty_t, empty_t>* neighbors[3];
    // Constructor
    TopologyCluster() : parent(), neighbors(){};
    // Helper functions
    int get_degree();
    bool contracts();
    bool contains_neighbor(TopologyCluster<empty_t, empty_t>* c);
    void insert_neighbor(TopologyCluster<empty_t, empty_t>* c);
    void remove_neighbor(TopologyCluster<empty_t, empty_t>* c);
    TopologyCluster<empty_t, empty_t>* get_root();
};

template<typename v_t, typename e_t>
class TopologyTree : ITernarizable {
using Cluster = TopologyCluster<v_t, e_t>;
public:
    // Topology tree interface
    TopologyTree(
        vertex_t n, QueryType q = CONNECTIVITY,
        std::function<v_t(v_t, v_t)> f_v = [](v_t x, v_t y) -> v_t {return x;},
        std::function<e_t(e_t, e_t)> f_e = [](e_t x, e_t y) -> e_t {return x;});
    TopologyTree(
        vertex_t n, QueryType q,
        std::function<v_t(v_t, v_t)> f_v, std::function<e_t(e_t, e_t)> f_e,
        v_t id_v, e_t id_e, v_t dval_v, e_t dval_e);
    ~TopologyTree();
    void link(vertex_t u, vertex_t v, e_t value);
    void link(vertex_t u, vertex_t v) { link(u,v,default_e); };
    void cut(vertex_t u, vertex_t v);
    void batch_link(Edge* links, int len);
    void batch_cut(Edge* cuts, int len);
    bool connected(vertex_t u, vertex_t v);
    v_t subtree_query(vertex_t v, vertex_t p = MAX_VERTEX_T);
    e_t path_query(vertex_t u, vertex_t v);
    //Interface methods overriden.
    short get_degree(vertex_t v) override {return leaves[v].get_degree();}
    std::pair<vertex_t, int> retrieve_v_to_del(vertex_t v) override {
      return std::pair(leaves[v].neighbors[0] - &leaves[0], 0); 
    }
    // Testing helpers
    bool is_valid();
    int get_height(vertex_t v);
    void print_tree();
    size_t space();
private:
    // Class data and parameters
    std::vector<Cluster> leaves;
    std::vector<std::vector<Cluster*>> root_clusters;
    std::vector<std::pair<std::pair<Cluster*, Cluster*>, bool>> contractions;
    QueryType query_type;
    std::function<v_t(v_t, v_t)> f_v;
    v_t identity_v;
    v_t default_v;
    std::function<e_t(e_t, e_t)> f_e;
    e_t identity_e;
    e_t default_e;
    // Helper functions
    void remove_ancestors(Cluster* c, int start_level = 0);
    void recluster_tree();
    void add_parent_adjacency(Cluster* c1, Cluster* c2, int i);
    void recompute_parent_value(Cluster* c1, Cluster* c2);
};

template<typename v_t>
class TopologyTree<v_t, empty_t> : ITernarizable {
using Cluster = TopologyCluster<v_t, empty_t>;
public:
    // Topology tree interface
    TopologyTree(
        vertex_t n, QueryType q = CONNECTIVITY,
        std::function<v_t(v_t, v_t)> f_v = [](v_t x, v_t y) -> v_t {return x;});
    TopologyTree(
        vertex_t n, QueryType q,
        std::function<v_t(v_t, v_t)> f_v, v_t id_v, v_t dval_v);
    ~TopologyTree();
    void link(vertex_t u, vertex_t v);
    void cut(vertex_t u, vertex_t v);
    bool connected(vertex_t u, vertex_t v);
    v_t subtree_query(vertex_t v, vertex_t p = MAX_VERTEX_T);
    //Interface methods overriden.
    short get_degree(vertex_t v) override {return leaves[v].get_degree();}
    std::pair<vertex_t, int> retrieve_v_to_del(vertex_t v) override {
      return std::pair(leaves[v].neighbors[0] - &leaves[0], 0); 
    }
    // Testing helpers
    bool is_valid();
    int get_height(vertex_t v);
    void print_tree();
    size_t space();
private:
    // Class data and parameters
    std::vector<Cluster> leaves;
    std::vector<std::vector<Cluster*>> root_clusters;
    std::vector<std::pair<std::pair<Cluster*, Cluster*>, bool>> contractions;
    QueryType query_type;
    std::function<v_t(v_t, v_t)> f_v;
    v_t identity_v;
    v_t default_v;
    // Helper functions
    void remove_ancestors(Cluster* c, int start_level = 0);
    void recluster_tree();
    void add_parent_adjacency(Cluster* c1, Cluster* c2, int i);
    void recompute_parent_value(Cluster* c1, Cluster* c2);
};

template<>
class TopologyTree<empty_t, empty_t> : ITernarizable {
using Cluster = TopologyCluster<empty_t, empty_t>;
public:
    // Topology tree interface
    TopologyTree(
        vertex_t n, QueryType q = CONNECTIVITY);
    ~TopologyTree();
    void link(vertex_t u, vertex_t v);
    void cut(vertex_t u, vertex_t v);
    bool connected(vertex_t u, vertex_t v);
    //Interface methods overriden.
    short get_degree(vertex_t v) override {return leaves[v].get_degree();}
    std::pair<vertex_t, int> retrieve_v_to_del(vertex_t v) override {
      return std::pair(leaves[v].neighbors[0] - &leaves[0], 0);
    }
    // Testing helpers
    bool is_valid();
    int get_height(vertex_t v);
    void print_tree();
    size_t space();
private:
    // Class data and parameters
    std::vector<Cluster> leaves;
    std::vector<std::vector<Cluster*>> root_clusters;
    std::vector<std::pair<std::pair<Cluster*, Cluster*>, bool>> contractions;
    QueryType query_type;
    // Helper functions
    void remove_ancestors(Cluster* c, int start_level = 0);
    void recluster_tree();
    void add_parent_adjacency(Cluster* c1, Cluster* c2, int i);
    void recompute_parent_value(Cluster* c1, Cluster* c2);
};

template<typename v_t, typename e_t>
TopologyTree<v_t, e_t>::TopologyTree(vertex_t n, QueryType q,
        std::function<v_t(v_t, v_t)> f_v, std::function<e_t(e_t, e_t)> f_e)
    : query_type(q), f_v(f_v), f_e(f_e) {
    leaves.resize(n);
    root_clusters.resize(max_tree_height(n));
    contractions.reserve(12);
}

template<typename v_t, typename e_t>
TopologyTree<v_t, e_t>::TopologyTree(vertex_t n, QueryType q,
        std::function<v_t(v_t, v_t)> f_v, std::function<e_t(e_t, e_t)> f_e,
        v_t id_v, e_t id_e, v_t dval_v, e_t dval_e)
    : query_type(q), f_v(f_v), f_e(f_e), identity_v(id_v), identity_e(id_e),
     default_v(dval_v), default_e(dval_e) {
    leaves.resize(n);
    root_clusters.resize(max_tree_height(n));
    contractions.reserve(12);
}

template<typename v_t, typename e_t>
TopologyTree<v_t, e_t>::~TopologyTree() {
    for (auto leaf : leaves) remove_ancestors(&leaf); // Clear all memory
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

template<typename v_t, typename e_t>
size_t TopologyTree<v_t, e_t>::space(){ 
    std::unordered_set<Cluster*> visited;
    size_t memory = sizeof(TopologyTree<v_t, e_t>);
    for(auto cluster : leaves){
        memory += sizeof(cluster);
        auto parent = cluster.parent;
        while(parent != nullptr && visited.count(parent) == 0){
            memory += sizeof(*parent);
            visited.insert(parent);
            parent = parent->parent;
        }
    }
    return memory;
}
/* Link vertex u and vertex v in the tree. Optionally include an
augmented value for the new edge (u,v). If no augmented value is
provided, the default value is 1. */
template<typename v_t, typename e_t>
void TopologyTree<v_t, e_t>::link(vertex_t u, vertex_t v, e_t value) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && !connected(u,v));
    START_TIMER(topology_remove_ancestor_timer);
    remove_ancestors(&leaves[u]);
    remove_ancestors(&leaves[v]);
    STOP_TIMER(topology_remove_ancestor_timer, topology_remove_ancestor_time);
    leaves[u].insert_neighbor(&leaves[v], value);
    leaves[v].insert_neighbor(&leaves[u], value);
    START_TIMER(topology_recluster_tree_timer);
    recluster_tree();
    STOP_TIMER(topology_recluster_tree_timer, topology_recluster_tree_time);
    // Collect tree height stats at the end of each update
    #ifdef COLLECT_HEIGHT_STATS
        max_height = std::max(max_height, get_height(u));
        max_height = std::max(max_height, get_height(v));
    #endif
}
template<typename v_t>
void TopologyTree<v_t, empty_t>::link(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && !connected(u,v));
    START_TIMER(topology_remove_ancestor_timer);
    remove_ancestors(&leaves[u]);
    remove_ancestors(&leaves[v]);
    STOP_TIMER(topology_remove_ancestor_timer, topology_remove_ancestor_time);
    leaves[u].insert_neighbor(&leaves[v]);
    leaves[v].insert_neighbor(&leaves[u]);
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
template<typename v_t, typename e_t>
void TopologyTree<v_t, e_t>::cut(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(leaves[u].contains_neighbor(&leaves[v]));
    START_TIMER(topology_remove_ancestor_timer);
    remove_ancestors(&leaves[u]);
    remove_ancestors(&leaves[v]);
    STOP_TIMER(topology_remove_ancestor_timer, topology_remove_ancestor_time);
    leaves[u].remove_neighbor(&leaves[v]);
    leaves[v].remove_neighbor(&leaves[u]);
    START_TIMER(topology_recluster_tree_timer);
    recluster_tree();
    STOP_TIMER(topology_recluster_tree_timer, topology_recluster_tree_time);
    // Collect tree height stats at the end of each update
    #ifdef COLLECT_HEIGHT_STATS
        max_height = std::max(max_height, get_height(u));
        max_height = std::max(max_height, get_height(v));
    #endif
}

template<typename v_t, typename e_t>
void TopologyTree<v_t, e_t>::remove_ancestors(Cluster* c, int start_level) {
    int level = start_level;
    for (auto neighbor : c->neighbors) {
        if (neighbor && neighbor->parent == c->parent) {
            neighbor->parent = nullptr; // Set sibling parent pointer to null
      root_clusters[level].push_back(
         neighbor);   // Keep track of parentless cluster
        }
    }
    auto curr = c->parent;
    c->parent = nullptr;
    root_clusters[level].push_back(c);
    while (curr) {
        auto prev = curr;
        curr = prev->parent;
        level++;
        for (auto neighbor : prev->neighbors) {
            if (neighbor && neighbor->parent == prev->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
        root_clusters[level].push_back(
           neighbor);   // Keep track of parentless cluster
            }
      if (neighbor)
        neighbor->remove_neighbor(prev);   // Remove prev from adjacency
        }
    auto position = std::find(root_clusters[level].begin(),
                              root_clusters[level].end(), prev);
    if (position != root_clusters[level].end())
      root_clusters[level].erase(position);
        delete prev; // Remove cluster prev
    }
}

template<typename v_t, typename e_t>
void TopologyTree<v_t, e_t>::recluster_tree() {
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
        for (auto cluster : root_clusters[level]) {
            if (cluster->get_degree() == 3) {
                // Combine deg 3 root clusters with deg 1 root  or non-root clusters
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor && neighbor->get_degree() == 1) {
                        auto parent = neighbor->parent;
                        bool new_parent = (parent == nullptr);
                        if (new_parent) { // If neighbor is a root cluster
                            parent = new Cluster();
                            root_clusters[level+1].push_back(parent);
                        }
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        recompute_parent_value(cluster, neighbor);
                        contractions.push_back({{cluster,neighbor},new_parent});
                        break;
                    }
                }
            } else if (cluster->get_degree() == 2 && !cluster->parent) {
                // Combine deg 2 root clusters with deg 1 or 2 root clusters
                for (int i = 0; i < 3; i++) {
                    auto neighbor = cluster->neighbors[i];
                    if (neighbor && !neighbor->parent && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
                        auto parent = new Cluster();
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        recompute_parent_value(cluster, neighbor);
                        root_clusters[level+1].push_back(parent);
                        contractions.push_back({{cluster,neighbor},true});
                        break;
                    }
                }
                // Combine deg 2 root clusters with deg 1 or 2 non-root clusters
                if (!cluster->parent)
                    for (int i = 0; i < 3; i++) {
                        auto neighbor = cluster->neighbors[i];
                        if (neighbor && neighbor->parent && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
                            if (neighbor->contracts())
                                continue;
                            auto parent = neighbor->parent;
                            if (!parent)
                                parent = new Cluster();
                            cluster->parent = parent;
                            neighbor->parent = parent;
                            recompute_parent_value(cluster, neighbor);
                            contractions.push_back({{cluster,neighbor},false});
                            break;
                        }
                    }
            } else if (cluster->get_degree() == 1 && !cluster->parent) {
                // Combine deg 1 root clusters with deg 1 root or non-root clusters
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor && neighbor->get_degree() == 1) {
                        auto parent = neighbor->parent;
                        bool new_parent = (parent == nullptr);
                        if (new_parent) { // If neighbor is a root cluster
                            parent = new Cluster();
                            root_clusters[level+1].push_back(parent);
                        }
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        recompute_parent_value(cluster, neighbor);
                        contractions.push_back({{cluster,neighbor},new_parent});
                        break;
                    }
                }
                // Combine deg 1 root clusters with deg 2 or 3 non-root clusters
                if (!cluster->parent)
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor && neighbor->parent && (neighbor->get_degree() == 2 || neighbor->get_degree() == 3)) {
                        if (neighbor->contracts())
                            continue;
                        auto parent = neighbor->parent;
                        if (!parent)
                            parent = new Cluster();
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        recompute_parent_value(cluster, neighbor);
                        contractions.push_back({{cluster,neighbor},false});
                        break;
                    }
                }
            }
            // Add remaining uncombined clusters to the next level
            if (!cluster->parent && cluster->get_degree() > 0) {
                auto parent = new Cluster();
                parent->value = cluster->value;
                cluster->parent = parent;
                recompute_parent_value(cluster, cluster);
                root_clusters[level+1].push_back(parent);
                contractions.push_back({{cluster,cluster},true});
            }
        }
        // Fill in the neighbor lists of the new clusters
        for (auto contraction : contractions) {
            auto c1 = contraction.first.first;
            auto c2 = contraction.first.second;
            auto parent = c1->parent;
            bool new_parent = contraction.second;
            for (int i = 0; i < 3; ++i)
                parent->neighbors[i] = nullptr;
            for (int i = 0; i < 3; ++i)
                if (c1->neighbors[i] && c1->neighbors[i] != c2)   // Don't add c2's parent (itself)
                    add_parent_adjacency(c1, c1->neighbors[i], i);
            for (int i = 0; i < 3; ++i)
                if (c2->neighbors[i] && c2->neighbors[i] != c1)   // Don't add c1's parent (itself)
                    add_parent_adjacency(c2, c2->neighbors[i], i);
            if (!new_parent)
                remove_ancestors(parent, level + 1);
        }
        // Clear the contents of this level
        root_clusters[level].clear();
        contractions.clear();
    }
}

template<typename v_t, typename e_t>
void TopologyTree<v_t, e_t>::add_parent_adjacency(Cluster* c1, Cluster* c2, int i) {
    assert(c1->neighbors[i] == c2);
    assert(c1->parent != c2->parent);
    c1->parent->insert_neighbor(c2->parent, c1->edge_values[i]);
    c2->parent->insert_neighbor(c1->parent, c1->edge_values[i]);
}

template<typename v_t>
void TopologyTree<v_t, empty_t>::add_parent_adjacency(Cluster* c1, Cluster* c2, int i) {
    assert(c1->neighbors[i] == c2);
    assert(c1->parent != c2->parent);
    c1->parent->insert_neighbor(c2->parent);
    c2->parent->insert_neighbor(c1->parent);
}

void TopologyTree<empty_t, empty_t>::add_parent_adjacency(Cluster* c1, Cluster* c2, int i) {
    assert(c1->neighbors[i] == c2);
    assert(c1->parent != c2->parent);
    c1->parent->insert_neighbor(c2->parent);
    c2->parent->insert_neighbor(c1->parent);
}

template<typename v_t, typename e_t>
void TopologyTree<v_t, e_t>::recompute_parent_value(Cluster* c1, Cluster* c2) {
    assert(c1->parent == c2->parent);
    auto parent = c1->parent;
    if (c1 == c2) {
        parent->value = c1->value;
        return;
    }
    if (query_type == SUBTREE) {
        parent->value = f_v(c1->value, c2->value);
    } else if (query_type == PATH && c1->get_degree() == 2 && c2->get_degree() == 2) {
        e_t edge_val;
        for (int i = 0; i < 3; i++)
            if (c1->neighbors[i] == c2)
                edge_val = c1->edge_values[i];
        parent->value = f_e(f_e(c1->value, c2->value), edge_val);
    }
}

template<typename v_t>
void TopologyTree<v_t, empty_t>::recompute_parent_value(Cluster* c1, Cluster* c2) {
    assert(c1->parent == c2->parent);
    auto parent = c1->parent;
    if (c1 == c2) {
        parent->value = c1->value;
        return;
    }
    if (query_type == SUBTREE) {
        parent->value = f_v(c1->value, c2->value);
    }
}

void TopologyTree<empty_t, empty_t>::recompute_parent_value(Cluster* c1, Cluster* c2) {
    assert(c1->parent == c2->parent);
    return;
}

template<typename v_t, typename e_t>
int TopologyCluster<v_t, e_t>::get_degree() {
    int deg = 0;
  for (auto neighbor : this->neighbors)
    if (neighbor)
      deg++;
    return deg;
}

// Helper function which returns whether this cluster combines with another
// cluster.
template<typename v_t, typename e_t>
bool TopologyCluster<v_t, e_t>::contracts() {
    bool contracts = false;
    for (auto neighbor : this->neighbors)
        if (neighbor && neighbor->parent == this->parent)
            contracts = true;
    return contracts;
}

template<typename v_t, typename e_t>
bool TopologyCluster<v_t, e_t>::contains_neighbor(TopologyCluster<v_t, e_t>* c) {
  for (int i = 0; i < 3; ++i)
    if (this->neighbors[i] == c)
      return true;
    return false;
}

template<typename v_t, typename e_t>
void TopologyCluster<v_t, e_t>::insert_neighbor(TopologyCluster<v_t, e_t>* c, e_t value) {
    if (this->contains_neighbor(c)) return;
    for (int i = 0; i < 3; ++i) {
        if (this->neighbors[i] == nullptr) {
            this->neighbors[i] = c;
            this->edge_values[i] = value;
            return;
        }
    }
    throw std::invalid_argument("No space to insert neighbor"); 
}

template<typename v_t>
void TopologyCluster<v_t, empty_t>::insert_neighbor(TopologyCluster<v_t, empty_t>* c) {
    if (this->contains_neighbor(c)) return;
    for (int i = 0; i < 3; ++i) {
        if (this->neighbors[i] == nullptr) {
            this->neighbors[i] = c;
            return;
        }
    }
    throw std::invalid_argument("No space to insert neighbor"); 
}

template<typename v_t, typename e_t>
void TopologyCluster<v_t, e_t>::remove_neighbor(TopologyCluster<v_t, e_t>* c) {
    for (int i = 0; i < 3; ++i) {
        if (this->neighbors[i] == c) {
            this->neighbors[i] = nullptr;
            return;
        }
    }
    throw std::invalid_argument("Neighbor to delete not found");
}

template<typename v_t, typename e_t>
TopologyCluster<v_t, e_t>* TopologyCluster<v_t, e_t>::get_root() {
    TopologyCluster<v_t, e_t>* curr = this;
    while (curr->parent) curr = curr->parent;
    return curr;
}

/* Return true if and only if there is a path from vertex u to
vertex v in the tree. */
template<typename v_t, typename e_t>
bool TopologyTree<v_t, e_t>::connected(vertex_t u, vertex_t v) {
    return leaves[u].get_root() == leaves[v].get_root();
}

/* Returns the value of the associative function f applied over
the augmented values for all the vertices in the subtree rooted
at v with respect to its parent p. If p = -1 (MAX_VERTEX_T) then
return the sum over the entire tree containing v. */
template<typename v_t, typename e_t>
v_t TopologyTree<v_t, e_t>::subtree_query(vertex_t v, vertex_t p) {
  assert(v >= 0 && v < leaves.size() && p >= 0 &&
         (p < leaves.size() || p == MAX_VERTEX_T));
  if (p == MAX_VERTEX_T)
    return leaves[v].get_root()->value;
    assert(leaves[v].contains_neighbor(&leaves[p]));
    // Get the total up until the LCA of v and p
    v_t total = f_v(identity_v, leaves[v].value);
    auto curr_v = &leaves[v];
    auto curr_p = &leaves[p];
    while (curr_v->parent != curr_p->parent) {
    for (auto neighbor : curr_v->neighbors)
      if (neighbor && neighbor->parent == curr_v->parent)
        total =
           f_v(total, neighbor->value);   // Only count vertices on the side of v
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
             f_v(total, neighbor->value);   // Count all remaining root clusters
            curr_v = curr_v->parent;
        }
        return total;
    }
  // If the cluster of v was deg 2 when it combined, only count the clusters on
  // the side of v
    Cluster* curr_u;
  for (auto neighbor : curr_v->neighbors)
    if (neighbor && neighbor != curr_p)
        curr_u = neighbor; // Find the neighbor of curr_v that is not curr_p
    while (curr_v->parent) {
        if (curr_u->parent == curr_v->parent) {
      total = f_v(
         total,
         curr_u->value);   // Count the remaining root clusters on the side of v
      if (curr_u->get_degree() == 1)
        return total;
            if (curr_u->get_degree() == 3) {
                curr_v = curr_v->parent;
                while (curr_v->parent) {
          for (auto neighbor : curr_v->neighbors)
            if (neighbor && neighbor->parent == curr_v->parent)
              total = f_v(total,
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
template<typename v_t, typename e_t>
e_t TopologyTree<v_t, e_t>::path_query(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && connected(u,v));
    // Compute the path on both sides for both vertices until they combine
    e_t path_u1, path_u2, path_v1, path_v2;
    path_u1 = path_u2 = path_v1 = path_v2 = identity_e;
    Cluster *bdry_u1, *bdry_u2, *bdry_v1, *bdry_v2;
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
        // NOTE(ATHARVA): Make this all into one function.
        for (int i = 0; i < 3; i++) {
            auto neighbor = curr_u->neighbors[i];
            if (neighbor && neighbor->parent == curr_u->parent) {
                if (curr_u->get_degree() == 2) {
                    if (curr_u->parent->get_degree() == 2) {
                        // Binary to Binary
                        if (neighbor == bdry_u1) {
                            path_u1 = f_e(path_u1, f_e(curr_u->edge_values[i], neighbor->value));
                            bdry_u2 = bdry_u2->parent;
                            for (int i = 0; i < 3; i++)
                            if (curr_u->parent->neighbors[i] && curr_u->parent->neighbors[i] != bdry_u2)
                                bdry_u1 = curr_u->parent->neighbors[i];
                        } else {
                            path_u2 = f_e(path_u2, f_e(curr_u->edge_values[i], neighbor->value));
                            bdry_u1 = bdry_u1->parent;
                            for (int i = 0; i < 3; i++)
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
                        // Unary to Binary and degree 3 case
                        if (curr_u->get_degree() != 3) path_u1 = path_u2 = f_e(path_u1, curr_u->edge_values[i]);
                        bdry_u1 = curr_u->parent->neighbors[0] ? curr_u->parent->neighbors[0] : curr_u->parent->neighbors[1];
                        bdry_u2 = curr_u->parent->neighbors[2] ? curr_u->parent->neighbors[2] : curr_u->parent->neighbors[1];
                    } else {
                        // Unary to Unary
                        path_u1 = f_e(path_u1, f_e(curr_u->edge_values[i], neighbor->value));
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
        for (int i = 0; i < 3; i++) {
            auto neighbor = curr_v->neighbors[i];
            if (neighbor && neighbor->parent == curr_v->parent) {
                if (curr_v->get_degree() == 2) {
                    if (curr_v->parent->get_degree() == 2) {
                        // Binary to Binary
                        if (neighbor == bdry_v1) {
                            path_v1 = f_e(path_v1, f_e(curr_v->edge_values[i], neighbor->value));
                            bdry_v2 = bdry_v2->parent;
                            for (int i = 0; i < 3; i++)
                            if (curr_v->parent->neighbors[i] && curr_v->parent->neighbors[i] != bdry_v2)
                                bdry_v1 = curr_v->parent->neighbors[i];
                        } else {
                            path_v2 = f_e(path_v2, f_e(curr_v->edge_values[i], neighbor->value));
                            bdry_v1 = bdry_v1->parent;
                            for (int i = 0; i < 3; i++)
                            if (curr_v->parent->neighbors[i] && curr_v->parent->neighbors[i] != bdry_v1)
                                bdry_v2 = curr_v->parent->neighbors[i];
                        }
                    } else {
                        // Binary to Unary
                        path_v1 = (neighbor == bdry_v1) ? path_v2 : path_v1;
                    }
                } else {
                    if (curr_v->parent->get_degree() == 2) {
                        // Unary to Binary
                        if (curr_v->get_degree() != 3) 
                            path_v1 = path_v2 = f_e(path_v1, curr_v->edge_values[i]);
                        bdry_v1 = curr_v->parent->neighbors[0] ? curr_v->parent->neighbors[0] : curr_v->parent->neighbors[1];
                        bdry_v2 = curr_v->parent->neighbors[2] ? curr_v->parent->neighbors[2] : curr_v->parent->neighbors[1];
                    } else {
                        // Unary to Unary
                        path_v1 = f_e(path_v1, f_e(curr_v->edge_values[i], neighbor->value));
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
    e_t total = identity_e;
    if (curr_u->get_degree() == 2)
        total = f_e(total, (curr_v == bdry_u1) ? path_u1 : path_u2);
    else
        total = f_e(total, path_u1);
    if (curr_v->get_degree() == 2)
        total = f_e(total, (curr_u == bdry_v1) ? path_v1 : path_v2);
    else
        total = f_e(total, path_v1);
    // Add the value of the last edge
    for (int i = 0; i < 3; i++)
        if (curr_u->neighbors[i] == curr_v)
            total = f_e(total, curr_u->edge_values[i]);
    return total;
}

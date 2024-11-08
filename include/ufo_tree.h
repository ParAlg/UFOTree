#pragma once
#include "types.h"
#include "util.h"
#include <absl/container/flat_hash_set.h>
#include <absl/container/flat_hash_map.h>

/* These constants determines the maximum size of array of nieghbors and
the vector of neighbors for each UFOCluster. Any additional neighbors will
be stored in the hash set for efficiency. */
#define UFO_ARRAY_MAX 3

// #define COLLECT_ROOT_CLUSTER_STATS
#ifdef COLLECT_ROOT_CLUSTER_STATS
    static std::map<int, int> root_clusters_histogram;
#endif

struct UFOClusterBase {
    // UFO cluster data
    UFOClusterBase* parent = nullptr;
    UFOClusterBase* neighbors[UFO_ARRAY_MAX];
    // Constructor
    UFOClusterBase() :  parent(), neighbors(){};
    // Helper functions
    bool contracts() {
        for (auto neighbor : this->neighbors) if (neighbor && neighbor->parent == this->parent) return true;
        return false;
    }
    UFOClusterBase* get_neighbor() {
        for (auto neighbor : neighbors)
            if (neighbor) return neighbor;
        std::cerr << "No neighbor found to return." << std::endl;
        std::abort();
    }
    UFOClusterBase* get_root() {
        UFOClusterBase* curr = this;
        while (curr->parent) curr = curr->parent;
        return curr;
    }
};

template<typename v_t, typename e_t>
struct UFOCluster : public UFOClusterBase {
    // template type values
    absl::flat_hash_map<UFOClusterBase*, e_t>* neighbors_set = nullptr;
    e_t edge_values[UFO_ARRAY_MAX]; // Only for path queries
    v_t value;                      // Stores subtree values of cluster path values
    // Constructor
    UFOCluster() : UFOClusterBase(), edge_values(), value(){};
    UFOCluster(v_t value) : UFOClusterBase(), edge_values(), value(value) {};
    // Helper functions
    int get_degree() {
        int deg = 0;
        for (auto neighbor : this->neighbors) if (neighbor) deg++;
        if (neighbors_set) [[unlikely]] deg += neighbors_set->size();
        return deg;
    }
    bool parent_high_fanout() {
        assert(parent);
        int parent_degree = static_cast<UFOCluster<v_t, e_t>*>(parent)->get_degree();
        if (get_degree() == 1) {
            auto neighbor = static_cast<UFOCluster<v_t, e_t>*>(get_neighbor());
            if (neighbor->parent == parent)
            if (neighbor->get_degree() - parent_degree > 2) return true;
        } else {
            if (get_degree() - parent_degree > 2) return true;
        }
        return false;
    }
    bool contains_neighbor(UFOClusterBase* c) {
        for (auto neighbor : neighbors) if (neighbor && neighbor == c) return true;
        if (neighbors_set && neighbors_set->find(c) != neighbors_set->end()) return true;
        return false;
    }
    void insert_neighbor(UFOClusterBase* c) {
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) if (this->neighbors[i] == c) return;
        assert(!contains_neighbor(c));
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) {
            if (this->neighbors[i] == nullptr) {
                this->neighbors[i] = c;
                return;
            }
        }
        if (!neighbors_set)
            neighbors_set = new absl::flat_hash_map<UFOClusterBase*, e_t>();
        std::pair<UFOClusterBase*, e_t> insert_pair;
        insert_pair.first = c;
        neighbors_set->insert(insert_pair);
    }
    void insert_neighbor_with_value(UFOClusterBase* c, e_t value) {
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) if (this->neighbors[i] == c) return;
        assert(!contains_neighbor(c));
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) {
            if (this->neighbors[i] == nullptr) {
                this->neighbors[i] = c;
                this->edge_values[i] = value;
                return;
            }
        }
        if (!neighbors_set)
            neighbors_set = new absl::flat_hash_map<UFOClusterBase*, e_t>();
        neighbors_set->insert({c,value});
    }
    void remove_neighbor(UFOClusterBase* c) {
        assert(contains_neighbor(c));
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) {
            if (this->neighbors[i] == c) {
                this->neighbors[i] = nullptr;
                //this->edge_values[i] = id_e;
                if (neighbors_set) [[unlikely]] { // Put an element from the set into the array
                    auto replacement = *neighbors_set->begin();
                    this->neighbors[i] = replacement.first;
                    this->edge_values[i] = replacement.second;
                    neighbors_set->erase(replacement.first);
                    if (neighbors_set->empty()) {
                        delete neighbors_set;
                        neighbors_set = nullptr;
                    }
                }
                return;
            }
        }
        neighbors_set->erase(c);
        if (neighbors_set->empty()) {
            delete neighbors_set;
            neighbors_set = nullptr;
        }
    }
    size_t calculate_size(){
        size_t memory = sizeof(UFOCluster<v_t, e_t>);
        if (neighbors_set)
            memory += neighbors_set->bucket_count() * sizeof(std::pair<UFOClusterBase*, e_t>);
        return memory;
    }
};
template<typename v_t>
struct UFOCluster<v_t, empty_t> : public UFOClusterBase {
    // template type values
    absl::flat_hash_set<UFOClusterBase*>* neighbors_set = nullptr;
    v_t value;                      // Stores subtree values of cluster path values
    // Constructor
    UFOCluster() : UFOClusterBase(), value(){};
    // Helper functions
    int get_degree() {
        int deg = 0;
        for (auto neighbor : this->neighbors) if (neighbor) deg++;
        if (neighbors_set) [[unlikely]] deg += neighbors_set->size();
        return deg;
    }
    bool parent_high_fanout() {
        assert(parent);
        int parent_degree = static_cast<UFOCluster<v_t, empty_t>*>(parent)->get_degree();
        if (get_degree() == 1) {
            auto neighbor = static_cast<UFOCluster<v_t, empty_t>*>(get_neighbor());
            if (neighbor->parent == parent)
            if (neighbor->get_degree() - parent_degree > 2) return true;
        } else {
            if (get_degree() - parent_degree > 2) return true;
        }
        return false;
    }
    bool contains_neighbor(UFOClusterBase* c) {
        for (auto neighbor : neighbors) if (neighbor && neighbor == c) return true;
        if (neighbors_set && neighbors_set->find(c) != neighbors_set->end()) return true;
        return false;
    }
    void insert_neighbor(UFOClusterBase* c) {
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) if (this->neighbors[i] == c) return;
        assert(!contains_neighbor(c));
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) {
            if (this->neighbors[i] == nullptr) {
                this->neighbors[i] = c;
                return;
            }
        }
        if (!neighbors_set)
            neighbors_set = new absl::flat_hash_set<UFOClusterBase*>();
        neighbors_set->insert(c);
    }
    void remove_neighbor(UFOClusterBase* c) {
        assert(contains_neighbor(c));
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) {
            if (this->neighbors[i] == c) {
                this->neighbors[i] = nullptr;
                if (neighbors_set) [[unlikely]] { // Put an element from the set into the array
                    auto replacement = *neighbors_set->begin();
                    this->neighbors[i] = replacement;
                    neighbors_set->erase(replacement);
                    if (neighbors_set->empty()) {
                        delete neighbors_set;
                        neighbors_set = nullptr;
                    }
                }
                return;
            }
        }
        neighbors_set->erase(c);
        if (neighbors_set->empty()) {
            delete neighbors_set;
            neighbors_set = nullptr;
        }
    }
    size_t calculate_size(){
        size_t memory = sizeof(UFOCluster<v_t, empty_t>);
        if (neighbors_set)
            memory += neighbors_set->bucket_count() * sizeof(UFOClusterBase*);
        return memory;
    }
};
template<>
struct UFOCluster<empty_t, empty_t> : public UFOClusterBase {
    // template type values
    absl::flat_hash_set<UFOClusterBase*>* neighbors_set = nullptr;
    // Constructor
    UFOCluster() : UFOClusterBase(){};
    // Helper functions
    int get_degree() {
        int deg = 0;
        for (auto neighbor : this->neighbors) if (neighbor) deg++;
        if (neighbors_set) [[unlikely]] deg += neighbors_set->size();
        return deg;
    }
    bool parent_high_fanout() {
        assert(parent);
        int parent_degree = static_cast<UFOCluster<empty_t, empty_t>*>(parent)->get_degree();
        if (get_degree() == 1) {
            auto neighbor = static_cast<UFOCluster<empty_t, empty_t>*>(get_neighbor());
            if (neighbor->parent == parent)
            if (neighbor->get_degree() - parent_degree > 2) return true;
        } else {
            if (get_degree() - parent_degree > 2) return true;
        }
        return false;
    }
    bool contains_neighbor(UFOClusterBase* c) {
        for (auto neighbor : neighbors) if (neighbor && neighbor == c) return true;
        if (neighbors_set && neighbors_set->find(c) != neighbors_set->end()) return true;
        return false;
    }
    void insert_neighbor(UFOClusterBase* c) {
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) if (this->neighbors[i] == c) return;
        assert(!contains_neighbor(c));
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) {
            if (this->neighbors[i] == nullptr) {
                this->neighbors[i] = c;
                return;
            }
        }
        if (!neighbors_set)
            neighbors_set = new absl::flat_hash_set<UFOClusterBase*>();
        neighbors_set->insert(c);
    }
    void remove_neighbor(UFOClusterBase* c) {
        assert(contains_neighbor(c));
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) {
            if (this->neighbors[i] == c) {
                this->neighbors[i] = nullptr;
                if (neighbors_set) [[unlikely]] { // Put an element from the set into the array
                    auto replacement = *neighbors_set->begin();
                    this->neighbors[i] = replacement;
                    neighbors_set->erase(replacement);
                    if (neighbors_set->empty()) {
                        delete neighbors_set;
                        neighbors_set = nullptr;
                    }
                }
                return;
            }
        }
        neighbors_set->erase(c);
        if (neighbors_set->empty()) {
            delete neighbors_set;
            neighbors_set = nullptr;
        }
    }
    size_t calculate_size(){
        size_t memory = sizeof(UFOCluster<empty_t, empty_t>);
        if (neighbors_set)
            memory += neighbors_set->bucket_count() * sizeof(UFOClusterBase*);
        return memory;
    }
};

template<typename v_t, typename e_t>
class UFOTree {
using Cluster = UFOCluster<v_t, e_t>;
public:
    // UFO tree interface
    UFOTree(
        vertex_t n, QueryType q = CONNECTIVITY,
        std::function<v_t(v_t, v_t)> f_v = [](v_t x, v_t y) -> v_t {return x;},
        std::function<e_t(e_t, e_t)> f_e = [](e_t x, e_t y) -> e_t {return x;});
    UFOTree(
        vertex_t n, QueryType q,
        std::function<v_t(v_t, v_t)> f_v, std::function<e_t(e_t, e_t)> f_e,
        v_t id_v, e_t id_e, v_t dval_v, e_t dval_e);
    UFOTree(int n, QueryType q, std::function<v_t(v_t, v_t)> f, v_t id, v_t d_val);
    ~UFOTree();
    void link(vertex_t u, vertex_t v);
    void link(vertex_t u, vertex_t v, e_t value);
    void cut(vertex_t u, vertex_t v);
    bool connected(vertex_t u, vertex_t v);
    e_t path_query(vertex_t u, vertex_t v);
    // Testing helpers
    size_t space();
    size_t count_nodes();
    size_t get_height();
    bool is_valid();
    void print_tree();
private:
    // Class data and parameters
    std::vector<Cluster> leaves;
    std::vector<std::vector<Cluster*>> root_clusters;
    std::vector<std::pair<std::pair<Cluster*,Cluster*>,bool>> contractions;
    std::vector<std::pair<Cluster*,vertex_t>> lower_deg[2]; // lower_deg helps to identify clusters who became low degree during a deletion update
    QueryType query_type;
    std::function<v_t(v_t, v_t)> f_v;
    v_t identity_v;
    v_t default_v;
    std::function<e_t(e_t, e_t)> f_e;
    e_t identity_e;
    e_t default_e;
    // We preallocate UFO clusters and store unused clusters in free_clusters
    std::vector<Cluster*> free_clusters;
    Cluster* allocate_cluster();
    void free_cluster(Cluster* c);
    // Helper functions
    void remove_ancestors(Cluster* c, int start_level = 0);
    void recluster_tree();
    bool is_high_degree_or_high_fanout(Cluster* cluster, Cluster* child, int level);
    void disconnect_siblings(Cluster* c, int level);
    void insert_adjacency(Cluster* u, Cluster* v);
    void insert_adjacency(Cluster* u, Cluster* v, e_t value);
    void remove_adjacency(Cluster* u, Cluster* v);
};

template<typename v_t, typename e_t>
UFOTree<v_t, e_t>::UFOTree(vertex_t n, QueryType q,
        std::function<v_t(v_t, v_t)> f_v, std::function<e_t(e_t, e_t)> f_e)
    : query_type(q), f_v(f_v), f_e(f_e) {
    leaves.resize(n);
    root_clusters.resize(max_tree_height(n));
    contractions.reserve(12);
    for (int i = 0; i < n; ++i)
        free_clusters.push_back(new Cluster());
}

template<typename v_t, typename e_t>
UFOTree<v_t, e_t>::UFOTree(vertex_t n, QueryType q,
        std::function<v_t(v_t, v_t)> f_v, std::function<e_t(e_t, e_t)> f_e,
        v_t id_v, e_t id_e, v_t dval_v, e_t dval_e)
    : query_type(q), f_v(f_v), f_e(f_e), identity_v(id_v), identity_e(id_e),
     default_v(dval_v), default_e(dval_e) {
    leaves.resize(n, default_v);
    root_clusters.resize(max_tree_height(n));
    contractions.reserve(12);
    for (int i = 0; i < n; ++i)
        free_clusters.push_back(new Cluster());
}

template<typename v_t, typename e_t>
UFOTree<v_t, e_t>::UFOTree(int n, QueryType q,
        std::function<v_t(v_t, v_t)> f, v_t id, v_t d_val)
    : query_type(q), f_v(f), identity_v(id), default_v(d_val) {
    if constexpr (std::is_same<v_t,e_t>::value) {
        f_e = f;
        identity_e = id;
        default_e = d_val;
    }
    leaves.resize(n, default_v);
    root_clusters.resize(max_tree_height(n));
    contractions.reserve(12);
    for (int i = 0; i < n; ++i)
        free_clusters.push_back(new Cluster());
}

template<typename v_t, typename e_t>
UFOTree<v_t, e_t>::~UFOTree() {
    // Clear all memory
    std::unordered_set<Cluster*> clusters;
    for (auto leaf : leaves) {
        auto curr = leaf.parent;
        while (curr) {
            clusters.insert(static_cast<Cluster*>(curr));
            curr = curr->parent;
        }
    }
    for (auto cluster : clusters) delete cluster;
    for (auto cluster : free_clusters) delete cluster;
    #ifdef COLLECT_ROOT_CLUSTER_STATS
    std::cout << "Number of root clusters: Frequency" << std::endl;
        for (auto entry : root_clusters_histogram)
            std::cout << entry.first << "\t" << entry.second << std::endl;
    #endif
}

template<typename v_t, typename e_t>
UFOCluster<v_t, e_t>* UFOTree<v_t, e_t>::allocate_cluster() {
    if (!free_clusters.empty()) {
        auto c = free_clusters.back();
        free_clusters.pop_back();
        return c;
    }
    return new Cluster();
}

template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::free_cluster(UFOCluster<v_t, e_t>* c) {
    c->parent = nullptr;
    for (int i = 0; i < UFO_ARRAY_MAX; ++i)
        c->neighbors[i] = nullptr;
    if (c->neighbors_set) [[unlikely]] delete c->neighbors_set;
    c->neighbors_set = nullptr;
    free_clusters.push_back(c);
}

template<typename v_t, typename e_t>
size_t UFOTree<v_t, e_t>::space() {
    std::unordered_set<Cluster*> visited;
    size_t memory = sizeof(UFOTree<v_t, e_t>);
    for(auto cluster : leaves){
        memory += cluster.calculate_size();
        auto parent = static_cast<Cluster*>(cluster.parent);
        while(parent != nullptr && visited.count(parent) == 0){
            memory += parent->calculate_size();
            visited.insert(parent);
            parent = static_cast<Cluster*>(parent->parent);
        }
    }
    return memory;
}

template<typename v_t, typename e_t>
size_t UFOTree<v_t, e_t>::count_nodes() {
    std::unordered_set<Cluster*> visited;
    size_t node_count = 0;
    for(auto cluster : leaves){
        node_count += 1;
        auto parent = static_cast<Cluster*>(cluster.parent);
        while(parent != nullptr && visited.count(parent) == 0){
            node_count += 1;
            visited.insert(parent);
            parent = static_cast<Cluster*>(parent->parent);
        }
    }
    return node_count;
}

template<typename v_t, typename e_t>
size_t UFOTree<v_t, e_t>::get_height() {
    size_t max_height = 0;
    for (vertex_t v = 0; v < leaves.size(); ++v) {
        size_t height = 0;
        Cluster* curr = &leaves[v];
        while (curr) {
            height++;
            curr = static_cast<Cluster*>(curr->parent);
        }
        max_height = std::max(max_height, height);
    }
    return max_height;
}

/* Link vertex u and vertex v in the tree. Optionally include an
augmented value for the new edge (u,v). If no augmented value is
provided, the default value is 1. */
template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::link(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && !connected(u,v));
    remove_ancestors(&leaves[u]);
    remove_ancestors(&leaves[v]);
    insert_adjacency(&leaves[u], &leaves[v]);
    recluster_tree();
}
template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::link(vertex_t u, vertex_t v, e_t value) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && !connected(u,v));
    remove_ancestors(&leaves[u]);
    remove_ancestors(&leaves[v]);
    insert_adjacency(&leaves[u], &leaves[v], value);
    recluster_tree();
}

/* Cut vertex u and vertex v in the tree. */
template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::cut(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(leaves[u].contains_neighbor(&leaves[v]));
    auto curr_u = &leaves[u];
    auto curr_v = &leaves[v];
    while (curr_u != curr_v) {
        lower_deg[0].push_back({curr_u, curr_u->get_degree()-1});
        lower_deg[1].push_back({curr_v, curr_v->get_degree()-1});
        curr_u = static_cast<Cluster*>(curr_u->parent);
        curr_v = static_cast<Cluster*>(curr_v->parent);
    }
    remove_ancestors(&leaves[u]);
    remove_ancestors(&leaves[v]);
    lower_deg[0].clear();
    lower_deg[1].clear();
    remove_adjacency(&leaves[u], &leaves[v]);
    recluster_tree();
}

/* Removes the ancestors of cluster c that are not high degree nor
high fan-out and add them to root_clusters. */
template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::remove_ancestors(Cluster* c, int start_level) {
    int level = start_level; // level is always the level of cluster prev, 0 being the leaves
    auto prev = c;
    auto curr = static_cast<Cluster*>(c->parent);
    bool del = false;
    while (curr) {
        // Different cases for if curr will or will not be deleted later
        if (!is_high_degree_or_high_fanout(curr, prev, level)) [[likely]] { // We will delete curr next round
            disconnect_siblings(prev, level);
            if (del) [[likely]] { // Possibly delete prev
                for (auto neighbor : prev->neighbors)
                    if (neighbor) static_cast<Cluster*>(neighbor)->remove_neighbor(prev); // Remove prev from adjacency
                auto position = std::find(root_clusters[level].begin(), root_clusters[level].end(), prev);
                if (position != root_clusters[level].end()) root_clusters[level].erase(position);
                free_cluster(prev);
            } else [[unlikely]] {
                prev->parent = nullptr;
                root_clusters[level].push_back(prev);
            }
            del = true;
        } else [[unlikely]] { // We will not delete curr next round
            if (del) [[likely]] { // Possibly delete prev
                for (auto neighbor : prev->neighbors)
                    if (neighbor) static_cast<Cluster*>(neighbor)->remove_neighbor(prev); // Remove prev from adjacency
                auto position = std::find(root_clusters[level].begin(), root_clusters[level].end(), prev);
                if (position != root_clusters[level].end()) root_clusters[level].erase(position);
                free_cluster(prev);
            } else [[unlikely]] if (prev->get_degree() <= 1) {
                prev->parent = nullptr;
                root_clusters[level].push_back(prev);
            }
            del = false;
        }
        // Update pointers
        prev = curr;
        curr = static_cast<Cluster*>(prev->parent);
        level++;
    }
    // DO LAST DELETIONS
    if (del) [[likely]] { // Possibly delete prev
        for (auto neighbor : prev->neighbors)
            if (neighbor) static_cast<Cluster*>(neighbor)->remove_neighbor(prev); // Remove prev from adjacency
        auto position = std::find(root_clusters[level].begin(), root_clusters[level].end(), prev);
        if (position != root_clusters[level].end()) root_clusters[level].erase(position);
        free_cluster(prev);
    } else [[unlikely]] root_clusters[level].push_back(prev);
}

template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::recluster_tree() {
    for (int level = 0; level < root_clusters.size(); level++) {
        if (root_clusters[level].empty()) [[unlikely]] continue;
        // Update root cluster stats if we are collecting them
        #ifdef COLLECT_ROOT_CLUSTER_STATS
            if (root_clusters_histogram.find(root_clusters[level].size()) == root_clusters_histogram.end())
                root_clusters_histogram[root_clusters[level].size()] = 1;
            else
                root_clusters_histogram[root_clusters[level].size()] += 1;
        #endif
        // Merge deg exactly 3 root clusters with all of its deg 1 neighbors
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 3) [[unlikely]] {
                auto parent = allocate_cluster();
                if constexpr (!std::is_same<e_t, empty_t>::value) {
                    parent->value = identity_v;
                }
                cluster->parent = parent;
                root_clusters[level+1].push_back(parent);
                bool first_contraction = true;
                for (auto neighborp : cluster->neighbors) {
                    auto neighbor = static_cast<Cluster*>(neighborp);
                    if (neighbor && neighbor->get_degree() == 1) [[unlikely]] {
                        auto old_parent = neighbor->parent;
                        neighbor->parent = cluster->parent;
                        int lev = level+1;
                        while (old_parent) {
                            auto temp = old_parent;
                            old_parent = old_parent->parent;
                            auto position = std::find(root_clusters[lev].begin(), root_clusters[lev].end(), temp);
                            if (position != root_clusters[lev].end()) root_clusters[lev].erase(position);
                            free_cluster(static_cast<Cluster*>(temp));
                            lev++;
                        }
                        if (first_contraction) contractions.push_back({{neighbor, cluster},true});
                        first_contraction = false;
                    }
                }
                if (first_contraction) contractions.push_back({{cluster, cluster}, true});
            }
        }
        // Always combine deg 1 root clusters with its neighboring cluster
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 1) [[unlikely]] {
                auto neighbor = static_cast<Cluster*>(cluster->get_neighbor());  // Deg 1 cluster only one neighbor
                if (neighbor->parent) {
                    if (neighbor->get_degree() == 2 && neighbor->contracts()) continue;
                    cluster->parent = neighbor->parent;
                    contractions.push_back({{cluster,neighbor},false});
                } else {
                    auto parent = allocate_cluster();
                    if constexpr (!std::is_same<e_t, empty_t>::value) {
                        parent->value = identity_v;
                    }
                    cluster->parent = parent;
                    neighbor->parent = parent;
                    root_clusters[level+1].push_back(parent);
                    contractions.push_back({{cluster,neighbor},true});
                }
                // For deg exactly 3 neighbor make sure all deg 1 neighbors combine with it
                if (neighbor->get_degree() == 3) [[unlikely]] {
                    for (auto entryp : neighbor->neighbors) {
                        auto entry = static_cast<Cluster*>(entryp);
                        if (entry && entry->get_degree() == 1 && entry->parent != neighbor->parent) [[unlikely]] {
                            auto old_parent = static_cast<Cluster*>(entry->parent);
                            entry->parent = neighbor->parent;
                            contractions.push_back({{entry, neighbor},false});
                            static_cast<Cluster*>(neighbor->parent)->remove_neighbor(old_parent);
                            auto position = std::find(root_clusters[level+1].begin(), root_clusters[level+1].end(), old_parent);
                            if (position != root_clusters[level+1].end()) root_clusters[level+1].erase(position);
                            if (old_parent) free_cluster(old_parent);
                        }
                    }
                }
            }
        }
        // Combine deg 2 root clusters with deg 2 root clusters
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 2) [[unlikely]] {
                for (auto neighborp : cluster->neighbors) {
                    auto neighbor = static_cast<Cluster*>(neighborp);
                    if (neighbor && !neighbor->parent && (neighbor->get_degree() == 2)) [[unlikely]] {
                        auto parent = allocate_cluster();
                        if constexpr (!std::is_same<e_t, empty_t>::value) {
                            parent->value = identity_v;
                        }
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        root_clusters[level+1].push_back(parent);
                        contractions.push_back({{cluster,neighbor},true});
                        break;
                    }
                }
            }
        }
        // Combine deg 2 root clusters with deg 1 or 2 non-root clusters
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 2) [[unlikely]] {
                for (auto neighborp : cluster->neighbors) {
                    auto neighbor = static_cast<Cluster*>(neighborp);
                    if (neighbor && neighbor->parent && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) [[unlikely]] {
                        if (neighbor->contracts()) continue;
                        cluster->parent = neighbor->parent;
                        contractions.push_back({{cluster,neighbor},false}); // The order here is important
                        break;
                    }
                }
            }
        }
        // Add remaining uncombined clusters to the next level
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() > 0) [[unlikely]] {
                auto parent = allocate_cluster();
                if constexpr (!std::is_same<v_t, empty_t>::value) {
                    parent->value = cluster->value;
                }
                cluster->parent = parent;
                root_clusters[level+1].push_back(parent);
                contractions.push_back({{cluster,cluster},true});
            }
        }
        // Fill in the neighbor lists of the new clusters
        for (auto contraction : contractions) {
            auto c1 = contraction.first.first;
            auto c2 = contraction.first.second;
            auto parent = static_cast<Cluster*>(c1->parent);
            bool new_parent = contraction.second;
            if (new_parent) {
                if constexpr (!std::is_same<e_t, empty_t>::value) {
                    for (int i = 0; i < UFO_ARRAY_MAX; i++) { parent->neighbors[i] = nullptr; parent->edge_values[i] = identity_e;}
                }
                for (int i = 0; i < UFO_ARRAY_MAX; i++) {
                    auto neighbor = static_cast<Cluster*>(c1->neighbors[i]);
                    if (neighbor && neighbor != c2 && parent != neighbor->parent) {
                        if constexpr (std::is_same<e_t, empty_t>::value) {
                            parent->insert_neighbor(neighbor->parent);
                            static_cast<Cluster*>(neighbor->parent)->insert_neighbor(parent);
                        } else {
                            parent->insert_neighbor_with_value(neighbor->parent, c1->edge_values[i]);
                            static_cast<Cluster*>(neighbor->parent)->insert_neighbor_with_value(parent, c1->edge_values[i]);
                        }
                    }
                }
                if (c1->neighbors_set) [[unlikely]]
                for (auto neighbor_pair : *c1->neighbors_set) {
                    Cluster* neighbor;
                    if constexpr (std::is_same<e_t, empty_t>::value) {
                        neighbor = static_cast<Cluster*>(neighbor_pair);
                    } else {
                        neighbor = static_cast<Cluster*>(neighbor_pair.first);
                    }
                    if (neighbor && neighbor != c2 && parent != neighbor->parent) {
                        if constexpr (std::is_same<e_t, empty_t>::value) {
                            parent->insert_neighbor(neighbor->parent);
                            static_cast<Cluster*>(neighbor->parent)->insert_neighbor(parent);
                        } else {
                            parent->insert_neighbor_with_value(neighbor->parent, neighbor_pair.second);
                            static_cast<Cluster*>(neighbor->parent)->insert_neighbor_with_value(parent, neighbor_pair.second);
                        }
                    }
                }
                if (c1 != c2) [[likely]] {
                  for (int i = 0; i < UFO_ARRAY_MAX; i++) {
                      auto neighbor = static_cast<Cluster*>(c2->neighbors[i]);
                      if (neighbor && neighbor != c1 && parent != neighbor->parent) {
                        if constexpr (std::is_same<e_t, empty_t>::value) {
                          parent->insert_neighbor(neighbor->parent);
                          static_cast<Cluster*>(neighbor->parent)->insert_neighbor(parent);
                        } else {
                          parent->insert_neighbor_with_value(neighbor->parent, c2->edge_values[i]);
                          static_cast<Cluster*>(neighbor->parent)->insert_neighbor_with_value(parent, c2->edge_values[i]);
                        }
                      }
                  }
                  if  (c2->neighbors_set) [[unlikely]]
                  for (auto neighbor_pair : *c2->neighbors_set){
                      Cluster* neighbor;
                      if constexpr (std::is_same<e_t, empty_t>::value) {
                        neighbor = static_cast<Cluster*>(neighbor_pair);
                      } else {
                        neighbor = static_cast<Cluster*>(neighbor_pair.first);
                      }
                      if (neighbor && neighbor != c1 && parent != neighbor->parent) {
                        if constexpr (std::is_same<e_t, empty_t>::value) {
                          parent->insert_neighbor(neighbor->parent);
                          static_cast<Cluster*>(neighbor->parent)->insert_neighbor(parent);
                        } else {
                          parent->insert_neighbor_with_value(neighbor->parent, neighbor_pair.second);
                          static_cast<Cluster*>(neighbor->parent)->insert_neighbor_with_value(parent, neighbor_pair.second);
                        }
                      }
                  }
                }
            } else {
                // We ordered contractions so c2 is the one that had a parent already
                if (c1->get_degree() == 2) {
                    for (int i = 0; i < UFO_ARRAY_MAX; i++){
                        auto neighbor = static_cast<Cluster*>(c1->neighbors[i]);
                        if (neighbor && neighbor != c2)
                        if constexpr (std::is_same<e_t, empty_t>::value) {
                            insert_adjacency(parent, static_cast<Cluster*>(neighbor->parent));
                        } else {
                            insert_adjacency(parent, static_cast<Cluster*>(neighbor->parent), c1->edge_values[i]);
                        }
                    }
                }
                remove_ancestors(parent, level+1);
            }
            if constexpr (!std::is_same<e_t, empty_t>::value) {
                if (c1 != c2 && c1->get_degree() == 2 && c2->get_degree() == 2) {
                    for (int i = 0; i < UFO_ARRAY_MAX; i++) {
                        auto neighbor = static_cast<Cluster*>(c1->neighbors[i]);
                        if (neighbor && neighbor == c2) {
                            parent->value = f_e(c1->value, f_e(c2->value, c1->edge_values[i]));
                        }
                    }
                }
            }
        }
        // Clear the contents of this level
        root_clusters[level].clear();
        contractions.clear();
    }
}

template<typename v_t, typename e_t>
bool UFOTree<v_t, e_t>::is_high_degree_or_high_fanout(Cluster* cluster, Cluster* child, int level) {
    int cluster_degree;
    if (level+1 < lower_deg[0].size() && lower_deg[0][level+1].first == cluster) [[unlikely]] cluster_degree = lower_deg[0][level+1].second;
    else if (level+1 < lower_deg[1].size() && lower_deg[1][level+1].first == cluster) [[unlikely]] cluster_degree = lower_deg[1][level+1].second;
    else [[likely]] cluster_degree = cluster->get_degree();
    if (cluster_degree > 2) [[unlikely]] return true;
    int child_degree;
    if (level < lower_deg[0].size() && lower_deg[0][level].first == child) [[unlikely]] child_degree = lower_deg[0][level].second;
    else if (level < lower_deg[1].size() && lower_deg[1][level].first == child) [[unlikely]] child_degree = lower_deg[1][level].second;
    else [[likely]] child_degree = child->get_degree();
    if (child->get_degree() == 1) {
        auto neighbor = static_cast<Cluster*>(child->get_neighbor());
        if (neighbor->parent == cluster && neighbor->get_degree() - cluster->get_degree() > 2) [[unlikely]] return true;
    } else if (child_degree - cluster_degree > 2) [[unlikely]] return true;
    return false;
}

template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::insert_adjacency(Cluster* u, Cluster* v) {
    auto curr_u = u;
    auto curr_v = v;
    while (curr_u && curr_v && curr_u != curr_v) {
        curr_u->insert_neighbor(curr_v);
        curr_v->insert_neighbor(curr_u);
        curr_u = static_cast<Cluster*>(curr_u->parent);
        curr_v = static_cast<Cluster*>(curr_v->parent);
    }
}
template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::insert_adjacency(Cluster* u, Cluster* v, e_t value) {
    auto curr_u = u;
    auto curr_v = v;
    while (curr_u && curr_v && curr_u != curr_v) {
        curr_u->insert_neighbor_with_value(curr_v, value);
        curr_v->insert_neighbor_with_value(curr_u, value);
        curr_u = static_cast<Cluster*>(curr_u->parent);
        curr_v = static_cast<Cluster*>(curr_v->parent);
    }
}

template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::remove_adjacency(Cluster* u, Cluster* v) {
    auto curr_u = u;
    auto curr_v = v;
    while (curr_u && curr_v && curr_u != curr_v) {
        curr_u->remove_neighbor(curr_v);
        curr_v->remove_neighbor(curr_u);
        curr_u = static_cast<Cluster*>(curr_u->parent);
        curr_v = static_cast<Cluster*>(curr_v->parent);
    }
}

/* Helper function which takes a cluster c and the level of that cluster. The function
should find every cluster that shares a parent with c, disconnect it from their parent
and add it as a root cluster to be processed. */
template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::disconnect_siblings(Cluster* c, int level) {
    if (c->get_degree() == 1) {
        auto center = static_cast<Cluster*>(c->get_neighbor());
        if (c->parent != center->parent) return;
        assert(center->get_degree() <= 5);
        if (center->parent == c->parent) {
            for (auto neighborp : center->neighbors) {
                auto neighbor = static_cast<Cluster*>(neighborp);
                if (neighbor && neighbor->parent == c->parent && neighbor != c) {
                    neighbor->parent = nullptr; // Set sibling parent pointer to null
                    root_clusters[level].push_back(neighbor); // Keep track of root clusters
                }
            }
            if (center->neighbors_set) [[unlikely]]
            for (auto neighbor_pair : *center->neighbors_set){
                Cluster* neighbor;
                if constexpr (std::is_same<e_t, empty_t>::value) {
                    neighbor = static_cast<Cluster*>(neighbor_pair);
                } else {
                    neighbor = static_cast<Cluster*>(neighbor_pair.first);
                }
                if (neighbor && neighbor->parent == c->parent && neighbor != c) {
                    neighbor->parent = nullptr; // Set sibling parent pointer to null
                    root_clusters[level].push_back(neighbor); // Keep track of root clusters
                }
            }
            center->parent = nullptr;
            root_clusters[level].push_back(center);
        }
    } else {
        assert(c->get_degree() <= 5);
        for (auto neighborp : c->neighbors) {
            auto neighbor = static_cast<Cluster*>(neighborp);
            if (neighbor && neighbor->parent == c->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].push_back(neighbor); // Keep track of root clusters
            }
        }
        if (c->neighbors_set) [[unlikely]]
        for (auto neighbor_pair : *c->neighbors_set){
            Cluster* neighbor;
            if constexpr (std::is_same<e_t, empty_t>::value) {
                neighbor = static_cast<Cluster*>(neighbor_pair);
            } else {
                neighbor = static_cast<Cluster*>(neighbor_pair.first);
            }
            if (neighbor && neighbor->parent == c->parent && neighbor != c) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].push_back(neighbor); // Keep track of root clusters
            }
        }
    }
}

/* Return true if and only if there is a path from vertex u to
vertex v in the tree. */
template<typename v_t, typename e_t>
bool UFOTree<v_t, e_t>::connected(vertex_t u, vertex_t v) {
    return leaves[u].get_root() == leaves[v].get_root();
}

template<typename v_t, typename e_t>
e_t UFOTree<v_t, e_t>::path_query(vertex_t u, vertex_t v){
  assert(u < leaves.size()-1 && v > 0 && u < v && connected(u, v)); 

    e_t path_u1, path_u2, path_v1, path_v2;
    path_u1 = path_u2 = path_v1 = path_v2 = identity_e;
    Cluster *bdry_u1, *bdry_u2, *bdry_v1, *bdry_v2;
    bdry_u1 = bdry_u2 = bdry_v1 = bdry_v2 = nullptr;
    if (leaves[u].get_degree() == 2) {
        bdry_u1 = static_cast<Cluster*>(leaves[u].neighbors[0] ? leaves[u].neighbors[0] : leaves[u].neighbors[1]);
        bdry_u2 = static_cast<Cluster*>(leaves[u].neighbors[2] ? leaves[u].neighbors[2] : leaves[u].neighbors[1]);
    }
    if (leaves[v].get_degree() == 2) {
        bdry_v1 = static_cast<Cluster*>(leaves[v].neighbors[0] ? leaves[v].neighbors[0] : leaves[v].neighbors[1]);
        bdry_v2 = static_cast<Cluster*>(leaves[v].neighbors[2] ? leaves[v].neighbors[2] : leaves[v].neighbors[1]);
    }
    auto curr_u = static_cast<Cluster*>(&leaves[u]);
    auto curr_v = static_cast<Cluster*>(&leaves[v]);
    while (curr_u->parent != curr_v->parent) { 
        // NOTE(ATHARVA): Make this all into one function.
        for (int i = 0; i < UFO_ARRAY_MAX; i++) {
            auto neighbor = static_cast<Cluster*>(curr_u->neighbors[i]);
            if (neighbor && neighbor->parent == curr_u->parent) {
                if (curr_u->get_degree() == 2) {
                    if (static_cast<Cluster*>(curr_u->parent)->get_degree() == 2) {
                        // Binary to Binary
                        if (neighbor == bdry_u1) {
                            path_u1 = f_e(path_u1, f_e(curr_u->edge_values[i], neighbor->value));
                            bdry_u2 = static_cast<Cluster*>(bdry_u2->parent);
                            for (int i = 0; i < UFO_ARRAY_MAX; i++)
                                if(curr_u->parent->neighbors[i] && static_cast<Cluster*>(curr_u->parent->neighbors[i]) != bdry_u2)
                                    bdry_u1 = static_cast<Cluster*>(curr_u->parent->neighbors[i]);
                        } else {
                            path_u2 = f_e(path_u2, f_e(curr_u->edge_values[i], neighbor->value));
                            bdry_u1 = static_cast<Cluster*>(bdry_u1->parent);
                            for (int i = 0; i < UFO_ARRAY_MAX; i++)
                                if(curr_u->parent->neighbors[i] && static_cast<Cluster*>(curr_u->parent->neighbors[i]) != bdry_u1)
                                    bdry_u2 = static_cast<Cluster*>(curr_u->parent->neighbors[i]);
                        }
                    } else {
                        // Binary to Unary
                        path_u1 = (neighbor == bdry_u1) ? path_u2 : path_u1;
                    }
                } else if(curr_u->get_degree() > 2){
                    if(static_cast<Cluster*>(curr_u->parent)->get_degree() == 2){ 
                      // Superunary to Superunary/Binary
                      bdry_u1 = static_cast<Cluster*>(curr_u->parent->neighbors[0] ? curr_u->parent->neighbors[0] : curr_u->parent->neighbors[1]);
                      bdry_u2 = static_cast<Cluster*>(curr_u->parent->neighbors[2] ? curr_u->parent->neighbors[2] : curr_u->parent->neighbors[1]);
                    } 
                }else {
                    if (static_cast<Cluster*>(curr_u->parent)->get_degree() == 2) {
                        // Unary to Binary
                        path_u1 = path_u2 = f_e(path_u1, curr_u->edge_values[i]);
                        bdry_u1 = static_cast<Cluster*>(curr_u->parent->neighbors[0] ? curr_u->parent->neighbors[0] : curr_u->parent->neighbors[1]);
                        bdry_u2 = static_cast<Cluster*>(curr_u->parent->neighbors[2] ? curr_u->parent->neighbors[2] : curr_u->parent->neighbors[1]);
                    } else {
                        // Unary to Unary and Unary to Superunary.
                        path_u1 = f_e(path_u1, f_e(curr_u->edge_values[i], neighbor->value));
                    }
                }
                break;
            }
        }
        if (!curr_u->contracts()) {
            if (bdry_u1) bdry_u1 = static_cast<Cluster*>(bdry_u1->parent);
            if (bdry_u2) bdry_u2 = static_cast<Cluster*>(bdry_u2->parent);
        }
        curr_u = static_cast<Cluster*>(curr_u->parent);
        for (int i = 0; i < UFO_ARRAY_MAX; i++) {
            auto neighbor = static_cast<Cluster*>(curr_v->neighbors[i]);
            if (neighbor && neighbor->parent == curr_v->parent) {
                if (curr_v->get_degree() == 2) {
                    if (static_cast<Cluster*>(curr_v->parent)->get_degree() == 2) {
                        // Binary to Binary
                        if (neighbor == bdry_v1) {
                            path_v1 = f_e(path_v1, f_e(curr_v->edge_values[i], neighbor->value));
                            bdry_v2 = static_cast<Cluster*>(bdry_v2->parent);
                            for (int i = 0; i < UFO_ARRAY_MAX; i++)
                                if(curr_v->parent->neighbors[i] && static_cast<Cluster*>(curr_v->parent->neighbors[i]) != bdry_v2)
                                    bdry_v1 = static_cast<Cluster*>(curr_v->parent->neighbors[i]);
                        } else {
                            path_v2 = f_e(path_v2, f_e(curr_v->edge_values[i], neighbor->value));
                            bdry_v1 = static_cast<Cluster*>(bdry_v1->parent);
                            for (int i = 0; i < UFO_ARRAY_MAX; i++)
                                if(curr_v->parent->neighbors[i] && static_cast<Cluster*>(curr_v->parent->neighbors[i]) != bdry_v1)
                                    bdry_v2 = static_cast<Cluster*>(curr_v->parent->neighbors[i]);
                        }
                    } else {
                        // Binary to Unary
                        path_v1 = (neighbor == bdry_v1) ? path_v2 : path_v1;
                    }
                } else if(curr_v->get_degree() > 2){
                    if(static_cast<Cluster*>(curr_v->parent)->get_degree() == 2){ 
                      // Superunary to Superunary/Binary
                      bdry_v1 = static_cast<Cluster*>(curr_v->parent->neighbors[0] ? curr_v->parent->neighbors[0] : curr_v->parent->neighbors[1]);
                      bdry_v2 = static_cast<Cluster*>(curr_v->parent->neighbors[2] ? curr_v->parent->neighbors[2] : curr_v->parent->neighbors[1]);
                    } 
                }else {
                    if (static_cast<Cluster*>(curr_v->parent)->get_degree() == 2) {
                        // Unary to Binary
                        if(curr_v->get_degree() != 3) 
                          path_v1 = path_v2 = f_e(path_v1, curr_v->edge_values[i]);
                        bdry_v1 = static_cast<Cluster*>(curr_v->parent->neighbors[0] ? curr_v->parent->neighbors[0] : curr_v->parent->neighbors[1]);
                        bdry_v2 = static_cast<Cluster*>(curr_v->parent->neighbors[2] ? curr_v->parent->neighbors[2] : curr_v->parent->neighbors[1]);
                    } else {
                        // Unary to Unary
                        path_v1 = f_e(path_v1, f_e(curr_v->edge_values[i], neighbor->value));
                    }
                }
                break;
            }
        }
        if (!curr_v->contracts()) {
            if (bdry_v1) bdry_v1 = static_cast<Cluster*>(bdry_v1->parent);
            if (bdry_v2) bdry_v2 = static_cast<Cluster*>(bdry_v2->parent);
        }
        curr_v = static_cast<Cluster*>(curr_v->parent);
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
    for (int i = 0; i < UFO_ARRAY_MAX; i++) if (static_cast<Cluster*>(curr_u->neighbors[i]) == curr_v)
        total = f_e(total, curr_u->edge_values[i]);
    return total;
}

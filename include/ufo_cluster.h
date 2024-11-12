#pragma once
#include "types.h"
#include "util.h"
#include <absl/container/flat_hash_set.h>
#include <absl/container/flat_hash_map.h>

/* These constants determines the maximum size of array of nieghbors and
the vector of neighbors for each UFOCluster. Any additional neighbors will
be stored in the hash set for efficiency. Minimum value is 3 for queries
function correctly. */
#define UFO_ARRAY_MAX 3

// #define COLLECT_ROOT_CLUSTER_STATS
#ifdef COLLECT_ROOT_CLUSTER_STATS
    static std::map<int, int> root_clusters_histogram;
#endif


template<typename v_t, typename e_t>
class UFOCluster {
using Cluster = UFOCluster<v_t, e_t>;
public:
    // UFO cluster data
    [[no_unique_address]] e_t edge_value1; // Note: [[no_unique_address]] fields must come first
    [[no_unique_address]] e_t edge_value2;
    [[no_unique_address]] e_t edge_value3;
    [[no_unique_address]] v_t value;
    Cluster* parent = nullptr;
    Cluster* neighbors[UFO_ARRAY_MAX];
    absl::flat_hash_map<Cluster*,e_t>* neighbors_set = nullptr;
    // Constructor
    UFOCluster() : parent(), neighbors(), edge_value1(), edge_value2(), edge_value3(), value() {};
    UFOCluster(v_t val) : parent(), neighbors(), edge_value1(), edge_value2(), edge_value3(), value(val) {};
    // Helper functions
    Cluster* get_neighbor();
    Cluster* get_root();
    bool contracts();
    int get_degree();
    bool parent_high_fanout();
    bool contains_neighbor(Cluster* c);
    void insert_neighbor(Cluster* c);
    void insert_neighbor_with_value(Cluster* c, e_t value);
    void remove_neighbor(Cluster* c);
    void set_edge_value(int index, e_t value);
    e_t get_edge_value(int index);
    size_t calculate_size();
};

template<typename v_t, typename e_t>
inline UFOCluster<v_t,e_t>* UFOCluster<v_t,e_t>::get_neighbor() {
    assert(UFO_ARRAY_MAX > 0 && neighbors[0]);
    return neighbors[0];
}

template<typename v_t, typename e_t>
inline UFOCluster<v_t,e_t>* UFOCluster<v_t,e_t>::get_root() {
    Cluster* curr = this;
    while (curr->parent) curr = curr->parent;
    return curr;
}

template<typename v_t, typename e_t>
inline bool UFOCluster<v_t,e_t>::contracts() {
    assert(get_degree() <= UFO_ARRAY_MAX);
    for (auto neighbor : this->neighbors) if (neighbor && neighbor->parent == this->parent) return true;
    return false;
}

template<typename v_t, typename e_t>
inline int UFOCluster<v_t,e_t>::get_degree() {
    int deg = 0;
    for (auto neighbor : this->neighbors) if (neighbor) deg++;
    if (neighbors_set) [[unlikely]] deg += neighbors_set->size();
    return deg;
}

template<typename v_t, typename e_t>
inline bool UFOCluster<v_t,e_t>::parent_high_fanout() {
    assert(parent);
    int parent_degree = parent->get_degree();
    if (get_degree() == 1) {
        auto neighbor = get_neighbor();
        if (neighbor->parent == parent)
        if (neighbor->get_degree() - parent_degree > 2) return true;
    } else {
        if (get_degree() - parent_degree > 2) return true;
    }
    return false;
}

template<typename v_t, typename e_t>
inline bool UFOCluster<v_t,e_t>::contains_neighbor(Cluster* c) {
    for (auto neighbor : neighbors) if (neighbor && neighbor == c) return true;
    if (neighbors_set && neighbors_set->find(c) != neighbors_set->end()) return true;
    return false;
}

template<typename v_t, typename e_t>
inline void UFOCluster<v_t,e_t>::insert_neighbor(Cluster* c) {
    assert(UFO_ARRAY_MAX >= 3); // Can we optimize this part out?
    for (int i = 0; i < UFO_ARRAY_MAX; ++i) if (this->neighbors[i] == c) return;
    assert(!contains_neighbor(c));
    for (int i = 0; i < UFO_ARRAY_MAX; ++i) {
        if (this->neighbors[i] == nullptr) [[likely]] {
            this->neighbors[i] = c;
            return;
        }
    }
    if (!neighbors_set)
        neighbors_set = new absl::flat_hash_map<Cluster*,e_t>();
    std::pair<Cluster*,e_t> insert_pair;
    insert_pair.first = c;
    neighbors_set->insert(insert_pair);
}

template<typename v_t, typename e_t>
inline void UFOCluster<v_t,e_t>::insert_neighbor_with_value(Cluster* c, e_t value) {
    if constexpr (!std::is_same<e_t, empty_t>::value) {
        assert(UFO_ARRAY_MAX >= 3);
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) if (this->neighbors[i] == c) return;
        assert(!contains_neighbor(c));
        for (int i = 0; i < UFO_ARRAY_MAX; ++i) {
            if (this->neighbors[i] == nullptr) [[likely]] {
                this->neighbors[i] = c;
                this->set_edge_value(i, value);
                return;
            }
        }
        if (!neighbors_set)
            neighbors_set = new absl::flat_hash_map<Cluster*, e_t>();
        neighbors_set->insert({c,value});
    }
}

template<typename v_t, typename e_t>
inline void UFOCluster<v_t,e_t>::remove_neighbor(Cluster* c) {
    assert(contains_neighbor(c));
    for (int i = 0; i < UFO_ARRAY_MAX; ++i) {
        if (this->neighbors[i] == c) {
            this->neighbors[i] = nullptr;
            if (neighbors_set) [[unlikely]] { // Put an element from the set into the array
                auto replacement = *neighbors_set->begin();
                this->neighbors[i] = replacement.first;
                if constexpr (!std::is_same<e_t,empty_t>::value)
                    this->set_edge_value(i, replacement.second);
                neighbors_set->erase(replacement.first);
                if (neighbors_set->empty()) {
                    delete neighbors_set;
                    neighbors_set = nullptr;
                }
            } else [[likely]] {
                for (int j = UFO_ARRAY_MAX-1; j > i; --j) {
                    if (neighbors[j]) [[unlikely]] {
                        neighbors[i] = neighbors[j];
                        neighbors[j] = nullptr;
                        if constexpr (!std::is_same<e_t,empty_t>::value)
                            set_edge_value(i, get_edge_value(j));
                        break;
                    }
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

template<typename v_t, typename e_t>
inline void UFOCluster<v_t,e_t>::set_edge_value(int index, e_t value) {
    e_t* address = &edge_value1 + index;
    *address = value;
}

template<typename v_t, typename e_t>
inline e_t UFOCluster<v_t,e_t>::get_edge_value(int index) {
    e_t* address = &edge_value1 + index;
    return *address;
}

template<typename v_t, typename e_t>
inline size_t UFOCluster<v_t,e_t>::calculate_size() {
    size_t memory = sizeof(UFOCluster<v_t, e_t>);
    if (neighbors_set) memory += neighbors_set->bucket_count() * sizeof(std::pair<Cluster*, e_t>);
    return memory;
}

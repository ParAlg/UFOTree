#pragma once
#include "types.h"
#include "util.h"
#include <absl/container/flat_hash_set.h>
#include <absl/container/flat_hash_map.h>

/* These constants determines the maximum size of array of neighbors  or
children for each UFOCluster. Any additional neighbors or children will
be stored in the hash set for efficiency. Minimum value is 3 for queries
function correctly. */
#define UFO_NEIGHBOR_MAX 3
#define UFO_CHILD_MAX 3

// #define COLLECT_ROOT_CLUSTER_STATS
#ifdef COLLECT_ROOT_CLUSTER_STATS
    static std::map<int, int> root_clusters_histogram;
#endif


template<typename v_t, typename e_t>
class UFOCluster {
using Cluster = UFOCluster<v_t, e_t>;
using NeighborSet = absl::flat_hash_map<Cluster*,e_t>;
using ChildSet = absl::flat_hash_set<Cluster*>;
public:
    // Query fields, note that the [[no_unique_address]] fields must be declared first
    [[no_unique_address]] e_t edge_value1;
    [[no_unique_address]] e_t edge_value2;
    [[no_unique_address]] e_t edge_value3;
    [[no_unique_address]] v_t value;
    // Parent pointer
    Cluster* parent = nullptr;
    /* We tag the last neighbor pointer in the array with information about the degree of the cluster.
    If it is 1, 2, or 3, that is the degree of the cluster. If it is 4, then the cluster has degree 4
    or higher and the last neighbor pointer is actually a pointer to the NeighborSet object containing
    the remaining neighbors of the cluster. */
    Cluster* neighbors[UFO_NEIGHBOR_MAX];
    /* We tag the last child pointer in the array with information about the fanout of the cluster. If
    it is 1 or 2, that is the fanout of the cluster. If it is 3, then the cluster has fanout 3 or higher
    and the last child pointer is actually a pointer to the ChildSet object containing the remaining
    children of the cluster. */
    Cluster* children[UFO_CHILD_MAX];
    int degree = 0;
    int fanout = 0;
    // Constructors
    UFOCluster() : neighbors(), children() {};
    UFOCluster(v_t val) : neighbors(), children(), value(val) {};
    // Neighbors
    bool has_neighbor_set();
    NeighborSet* get_neighbor_set();
    void insert_neighbor(Cluster* c);
    void insert_neighbor_with_value(Cluster* c, e_t value);
    void remove_neighbor(Cluster* c);
    bool contains_neighbor(Cluster* c);
    // Children
    bool has_child_set();
    ChildSet* get_child_set();
    void insert_child(Cluster* c);
    void remove_child(Cluster* c);
    // Additional helper functions
    Cluster* get_root();
    bool contracts();
    int get_degree();
    void set_edge_value(int index, e_t value);
    e_t get_edge_value(int index);
    size_t calculate_size();
};

template<typename v_t, typename e_t>
bool UFOCluster<v_t,e_t>::has_neighbor_set() {
    int tag = GET_TAG(neighbors[UFO_NEIGHBOR_MAX-1]);
    if (tag <= UFO_NEIGHBOR_MAX) [[likely]] return false;
    return true;
}

template<typename v_t, typename e_t>
absl::flat_hash_map<UFOCluster<v_t, e_t>*,e_t>* UFOCluster<v_t,e_t>::get_neighbor_set() {
    return (NeighborSet*) UNTAG(neighbors[UFO_NEIGHBOR_MAX-1]);
}

template<typename v_t, typename e_t>
bool UFOCluster<v_t,e_t>::contains_neighbor(Cluster* c) {
    for (auto neighbor : neighbors) if (UNTAG(neighbor) == c) return true;
    if (has_neighbor_set() && get_neighbor_set()->find(c) != get_neighbor_set()->end()) return true;
    return false;
}

template<typename v_t, typename e_t>
void UFOCluster<v_t,e_t>::insert_neighbor(Cluster* c) {
    assert(!contains_neighbor(c));
    // degree++;
    for (int i = 0; i < UFO_NEIGHBOR_MAX; ++i) {
        if (UNTAG(neighbors[i]) == nullptr) [[likely]] {
            int tag = GET_TAG(neighbors[UFO_NEIGHBOR_MAX-1]);
            neighbors[i] = c;
            neighbors[UFO_NEIGHBOR_MAX-1] = TAG(UNTAG(neighbors[UFO_NEIGHBOR_MAX-1]), tag+1);
            return;
        }
    }
    if (!has_neighbor_set()) {
        auto neighbor_set = new NeighborSet();
        std::pair<Cluster*,e_t> insert_pair;
        insert_pair.first = UNTAG(neighbors[UFO_NEIGHBOR_MAX-1]);
        neighbor_set->insert(insert_pair);
        neighbors[UFO_NEIGHBOR_MAX-1] = TAG(neighbor_set, UFO_NEIGHBOR_MAX+1);
    }
    std::pair<Cluster*,e_t> insert_pair;
    insert_pair.first = c;
    get_neighbor_set()->insert(insert_pair);
}

template<typename v_t, typename e_t>
void UFOCluster<v_t,e_t>::insert_neighbor_with_value(Cluster* c, e_t value) {
    if constexpr (!std::is_same<e_t, empty_t>::value) {
        assert(!contains_neighbor(c));
        // degree++;
        for (int i = 0; i < UFO_NEIGHBOR_MAX; ++i) {
            if (UNTAG(neighbors[i]) == nullptr) [[likely]] {
                int tag = GET_TAG(neighbors[UFO_NEIGHBOR_MAX-1]);
                neighbors[i] = c;
                set_edge_value(i, value);
                neighbors[UFO_NEIGHBOR_MAX-1] = TAG(UNTAG(neighbors[UFO_NEIGHBOR_MAX-1]), tag+1);
                return;
            }
        }
        if (!has_neighbor_set()) {
            auto neighbor_set = new NeighborSet();
            neighbor_set->insert({UNTAG(neighbors[UFO_NEIGHBOR_MAX-1]), get_edge_value(UFO_NEIGHBOR_MAX-1)});
            neighbors[UFO_NEIGHBOR_MAX-1] = TAG(neighbor_set, UFO_NEIGHBOR_MAX+1);
        }
        get_neighbor_set()->insert({c,value});
    }
}

template<typename v_t, typename e_t>
void UFOCluster<v_t,e_t>::remove_neighbor(Cluster* c) {
    assert(contains_neighbor(c));
    // degree--;
    for (int i = 0; i < UFO_NEIGHBOR_MAX; ++i) {
        if (UNTAG(neighbors[i]) == c) {
            neighbors[i] = TAG(nullptr, GET_TAG(neighbors[i]));
            if (has_neighbor_set()) [[unlikely]] { // Put an element from the set into the array
                auto neighbor_set = get_neighbor_set();
                auto replacement = *neighbor_set->begin();
                neighbors[i] = replacement.first;
                if constexpr (!std::is_same<e_t,empty_t>::value)
                    set_edge_value(i, replacement.second);
                neighbor_set->erase(replacement.first);
                if (neighbor_set->size() == 1) {
                    auto replacement = *neighbor_set->begin();
                    neighbors[UFO_NEIGHBOR_MAX-1] = TAG(replacement.first, UFO_NEIGHBOR_MAX);
                    if constexpr (!std::is_same<e_t,empty_t>::value)
                        set_edge_value(UFO_NEIGHBOR_MAX-1, replacement.second);
                    delete neighbor_set;
                }
            } else [[likely]] {
                for (int j = UFO_NEIGHBOR_MAX-1; j > i; --j) {
                    if (UNTAG(neighbors[j])) [[unlikely]] {
                        neighbors[i] = UNTAG(neighbors[j]);
                        neighbors[j] = TAG(nullptr, GET_TAG(neighbors[j]));
                        if constexpr (!std::is_same<e_t,empty_t>::value)
                            set_edge_value(i, get_edge_value(j));
                        break;
                    }
                }
                neighbors[UFO_NEIGHBOR_MAX-1] = TAG(UNTAG(neighbors[UFO_NEIGHBOR_MAX-1]), GET_TAG(neighbors[UFO_NEIGHBOR_MAX-1])-1);
            }
            return;
        }
    }
    auto neighbor_set = get_neighbor_set();
    neighbor_set->erase(c);
    if (neighbor_set->size() == 1) {
        auto replacement = *neighbor_set->begin();
        neighbors[UFO_NEIGHBOR_MAX-1] = TAG(replacement.first, UFO_NEIGHBOR_MAX);
        if constexpr (!std::is_same<e_t,empty_t>::value)
            set_edge_value(UFO_NEIGHBOR_MAX-1, replacement.second);
        delete neighbor_set;
    }
}

#define FOR_ALL_NEIGHBORS(C,F) {                            \
    if (!C->has_neighbor_set()) [[likely]] {                \
        for (int i = 0; i < UFO_NEIGHBOR_MAX; ++i) {        \
            Cluster* neighbor = UNTAG(C->neighbors[i]);     \
            if (neighbor) F(neighbor,C->get_edge_value(i)); \
        }                                                   \
    } else [[unlikely]] {                                   \
        for (int i = 0; i < UFO_NEIGHBOR_MAX-1; ++i) {      \
            Cluster* neighbor = C->neighbors[i];            \
            if (neighbor) F(neighbor,C->get_edge_value(i)); \
        }                                                   \
        for (auto neighbor_pair : *C->get_neighbor_set()) { \
            Cluster* neighbor = neighbor_pair.first;        \
            F(neighbor,neighbor_pair.second);               \
        }                                                   \
    }                                                       \
}

template<typename v_t, typename e_t>
bool UFOCluster<v_t,e_t>::has_child_set() {
    int tag = GET_TAG(children[UFO_CHILD_MAX-1]);
    if (tag <= UFO_CHILD_MAX) [[likely]] return false;
    return true;
}

template<typename v_t, typename e_t>
absl::flat_hash_set<UFOCluster<v_t, e_t>*>* UFOCluster<v_t,e_t>::get_child_set() {
    return (ChildSet*) UNTAG(children[UFO_CHILD_MAX-1]);
}

template<typename v_t, typename e_t>
void UFOCluster<v_t,e_t>::insert_child(Cluster* c) {
    fanout++;
    for (int i = 0; i < UFO_CHILD_MAX; ++i) {
        if (UNTAG(children[i]) == nullptr) [[likely]] {
            int tag = GET_TAG(children[UFO_CHILD_MAX-1]);
            children[i] = c;
            children[UFO_CHILD_MAX-1] = TAG(UNTAG(children[UFO_CHILD_MAX-1]), tag+1);
            return;
        }
    }
    if (!has_child_set()) {
        auto child_set = new ChildSet();
        child_set->insert(UNTAG(children[UFO_CHILD_MAX-1]));
        children[UFO_CHILD_MAX-1] = TAG(child_set, UFO_CHILD_MAX+1);
    }
    get_child_set()->insert(c);
}

template<typename v_t, typename e_t>
void UFOCluster<v_t,e_t>::remove_child(Cluster* c) {
    fanout--;
    for (int i = 0; i < UFO_CHILD_MAX; ++i) {
        if (UNTAG(children[i]) == c) {
            children[i] = TAG(nullptr, GET_TAG(children[i]));
            if (has_child_set()) [[unlikely]] { // Put an element from the set into the array
                auto child_set = get_child_set();
                auto replacement = *child_set->begin();
                children[i] = replacement;
                child_set->erase(replacement);
                if (child_set->size() == 1) {
                    children[UFO_CHILD_MAX-1] = TAG(*child_set->begin(), UFO_CHILD_MAX);
                    delete child_set;
                }
            } else [[likely]] {
                for (int j = UFO_CHILD_MAX-1; j > i; --j) {
                    if (UNTAG(children[j])) [[unlikely]] {
                        children[i] = UNTAG(children[j]);
                        children[j] = TAG(nullptr, GET_TAG(children[j]));
                        break;
                    }
                }
                children[UFO_CHILD_MAX-1] = TAG(UNTAG(children[UFO_CHILD_MAX-1]), GET_TAG(children[UFO_CHILD_MAX-1])-1);
            }
            return;
        }
    }
    auto child_set = get_child_set();
    child_set->erase(c);
    if (child_set->size() == 1) {
        children[UFO_CHILD_MAX-1] = TAG(*child_set->begin(), UFO_CHILD_MAX);
        delete child_set;
    }
}

#define FOR_ALL_CHILDREN(C,F) {                             \
    if (!C->has_child_set()) [[likely]] {                   \
        for (int i = 0; i < UFO_CHILD_MAX; ++i) {           \
            Cluster* child = UNTAG(C->children[i]);         \
            if (child) F(child);                            \
        }                                                   \
    } else [[unlikely]] {                                   \
        for (int i = 0; i < UFO_CHILD_MAX-1; ++i) {         \
            Cluster* child = C->children[i];                \
            if (child) F(child);                            \
        }                                                   \
        for (auto child : *C->get_child_set()) {            \
            F(child);                                       \
        }                                                   \
    }                                                       \
}

template<typename v_t, typename e_t>
UFOCluster<v_t,e_t>* UFOCluster<v_t,e_t>::get_root() {
    Cluster* curr = this;
    while (curr->parent) curr = curr->parent;
    return curr;
}

template<typename v_t, typename e_t>
bool UFOCluster<v_t,e_t>::contracts() {
    assert(get_degree() <= UFO_NEIGHBOR_MAX);
    for (auto neighborp : neighbors) {
        auto neighbor = UNTAG(neighborp);
        if (neighbor && neighbor->parent == parent) return true;
    }
    return false;
}

template<typename v_t, typename e_t>
int UFOCluster<v_t,e_t>::get_degree() {
    int tag = GET_TAG(neighbors[UFO_NEIGHBOR_MAX-1]);
    if (tag <= 3) [[likely]] return tag;
    return 2 + get_neighbor_set()->size();
}

template<typename v_t, typename e_t>
void UFOCluster<v_t,e_t>::set_edge_value(int index, e_t value) {
    e_t* address = &edge_value1 + index;
    *address = value;
}

template<typename v_t, typename e_t>
e_t UFOCluster<v_t,e_t>::get_edge_value(int index) {
    e_t* address = &edge_value1 + index;
    return *address;
}

template<typename v_t, typename e_t>
size_t UFOCluster<v_t,e_t>::calculate_size() {
    size_t memory = sizeof(UFOCluster<v_t, e_t>);
    if (has_neighbor_set()) memory += get_neighbor_set()->bucket_count() * sizeof(std::pair<Cluster*, e_t>);
    return memory;
}

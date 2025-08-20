#pragma once
#include <mutex>
#include "types.h"
#include "util.h"
#include "pam/pam.h"


namespace dgbs {

template <typename aug_t>
struct ParallelUFOCluster;

template <typename aug_t>
struct UFONeighbor {
    using key_t = ParallelUFOCluster<aug_t>*;
    static inline bool comp(key_t a, key_t b) { return a < b;}
};


template <typename aug_t>
struct ParallelUFOCluster {
    using Cluster = ParallelUFOCluster;
    using ufo_pam_set = pam_set<UFONeighbor<aug_t>>;

    // Cluster data
    ufo_pam_set neighbors;
    Cluster* parent;

    // Extra fields for parallel updates
    Cluster* partner;
    int degree;

    // Constructor
    ParallelUFOCluster() : neighbors(), parent(nullptr), partner(nullptr), degree(0) {};

    // Functions
    std::mutex mtx;
    void insert_neighbor(Cluster* c);
    void delete_neighbor(Cluster* c);
    void insert_neighbors(parlay::sequence<Cluster*>& cs);
    void delete_neighbors(parlay::sequence<Cluster*>& cs);
    void insert_neighbors_sorted(parlay::sequence<Cluster*>& cs);
    void delete_neighbors_sorted(parlay::sequence<Cluster*>& cs);

    int get_degree();
    bool contracts();
    bool atomic_contracts();
    bool contains_neighbor(Cluster* c);
    Cluster* get_neighbor();
    Cluster* get_other_neighbor(Cluster* c);
    Cluster* get_root();


    template <class F>
    void for_all_neighbors(const F& f) {
        ufo_pam_set::foreach_index(neighbors, [&] (auto cluster, size_t i) {
            f(cluster);
        });
    }

    parlay::sequence<std::pair<ParallelUFOCluster*, ParallelUFOCluster*>> get_filtered_in_edges() {
        return parlay::map_maybe(ufo_pam_set::entries(neighbors), [&] (auto neighbor) -> std::optional<std::pair<Cluster*, Cluster*>> {
            if (partner && !neighbor->partner) return std::make_pair(neighbor, this);
            return std::nullopt;
        });
    }

    void print_neighbors();
};

template <typename aug_t>
void ParallelUFOCluster<aug_t>::insert_neighbor(Cluster* c) {
    mtx.lock();
    neighbors = ufo_pam_set::insert(std::move(neighbors), c);
    degree = get_degree();
    mtx.unlock();
}

template <typename aug_t>
void ParallelUFOCluster<aug_t>::delete_neighbor(Cluster* c) {
    mtx.lock();
    neighbors = ufo_pam_set::remove(std::move(neighbors), c);
    degree = get_degree();
    mtx.unlock();
}

template <typename aug_t>
void ParallelUFOCluster<aug_t>::insert_neighbors(parlay::sequence<Cluster*>& cs) {
    neighbors = ufo_pam_set::multi_insert(std::move(neighbors), cs);
    degree = get_degree();
}

template <typename aug_t>
void ParallelUFOCluster<aug_t>::delete_neighbors(parlay::sequence<Cluster*>& cs) {
    neighbors = ufo_pam_set::multi_delete(std::move(neighbors), cs);
    degree = get_degree();
}

template <typename aug_t>
void ParallelUFOCluster<aug_t>::insert_neighbors_sorted(parlay::sequence<Cluster*>& cs) {
    neighbors = ufo_pam_set::multi_insert_sorted(std::move(neighbors), cs);
    degree = get_degree();
}

template <typename aug_t>
void ParallelUFOCluster<aug_t>::delete_neighbors_sorted(parlay::sequence<Cluster*>& cs) {
    neighbors = ufo_pam_set::multi_delete_sorted(std::move(neighbors), cs);
    degree = get_degree();
}

template <typename aug_t>
int ParallelUFOCluster<aug_t>::get_degree() {
    return neighbors.size();
}

template <typename aug_t>
bool ParallelUFOCluster<aug_t>::contracts() {
    if (!parent) return false;
    std::atomic<bool> contracts = false;
    for_all_neighbors([&] (auto neighbor) {
        if (neighbor->parent == parent)
            contracts = true;
    });
    return contracts;
}

template <typename aug_t>
bool ParallelUFOCluster<aug_t>::atomic_contracts() {
    auto parent = AtomicLoad(&this->parent);
    if (!parent) return false;
    std::atomic<bool> contracts = false;
    for_all_neighbors([&] (auto neighbor) {
        if (AtomicLoad(&neighbor->parent) == parent)
            contracts = true;
    });
    return contracts;
}

template<typename aug_t>
bool ParallelUFOCluster<aug_t>::contains_neighbor(Cluster* c) {
    return neighbors.contains(c);
}

template <typename aug_t>
ParallelUFOCluster<aug_t>* ParallelUFOCluster<aug_t>::get_neighbor() {
    if (neighbors.size() < 1) return nullptr;
    return *neighbors.select(0);
}

template <typename aug_t>
ParallelUFOCluster<aug_t>* ParallelUFOCluster<aug_t>::get_other_neighbor(Cluster* c) {
    if (neighbors.size() < 2) return nullptr;
    Cluster* neighbor1 = *neighbors.select(0);
    Cluster* neighbor2 = *neighbors.select(1);
    if (neighbor1 == c) return neighbor2;
    return neighbor1;
}

template <typename aug_t>
ParallelUFOCluster<aug_t>* ParallelUFOCluster<aug_t>::get_root() {
    Cluster* curr = this;
    while (curr->parent) curr = curr->parent;
    return curr;
}

}

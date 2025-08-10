#pragma once
#include "types.h"
#include "util.h"
#include <mutex>


namespace dgbs {

template <typename aug_t>
struct ParallelUFOCluster {
    // Topology cluster data
    std::unordered_set<ParallelUFOCluster*> neighbors;
    ParallelUFOCluster* parent;
    ParallelUFOCluster* partner;

    // Constructor
    ParallelUFOCluster() : parent(nullptr), partner(nullptr), neighbors() {};

    // Functions
    std::mutex mtx;
    void insert_neighbor(ParallelUFOCluster* c);
    void delete_neighbor(ParallelUFOCluster* c);
    void insert_neighbors(parlay::sequence<ParallelUFOCluster*>& cs);
    void delete_neighbors(parlay::sequence<ParallelUFOCluster*>& cs);

    bool contains_neighbor(ParallelUFOCluster* c);

    int get_degree();
    bool contracts();
    ParallelUFOCluster* get_neighbor();
    ParallelUFOCluster* get_other_neighbor(ParallelUFOCluster* c);
    ParallelUFOCluster* get_root();
    
    bool atomic_contracts();

    void print_neighbors();
};

template <typename aug_t>
void ParallelUFOCluster<aug_t>::insert_neighbor(ParallelUFOCluster* c) {
    mtx.lock();
    neighbors.insert(c);
    mtx.unlock();
}

template <typename aug_t>
void ParallelUFOCluster<aug_t>::delete_neighbor(ParallelUFOCluster* c) {
    mtx.lock();
    neighbors.erase(c);
    mtx.unlock();
}

template <typename aug_t>
void ParallelUFOCluster<aug_t>::insert_neighbors(parlay::sequence<ParallelUFOCluster*>& cs) {
    mtx.lock();
    for (auto c : cs) neighbors.insert(c);
    mtx.unlock();
}

template <typename aug_t>
void ParallelUFOCluster<aug_t>::delete_neighbors(parlay::sequence<ParallelUFOCluster*>& cs) {
    mtx.lock();
    for (auto c : cs) neighbors.erase(c);
    mtx.unlock();
}


template <typename aug_t>
bool ParallelUFOCluster<aug_t>::contains_neighbor(ParallelUFOCluster* c) {
    return neighbors.contains(c);
}

template <typename aug_t>
int ParallelUFOCluster<aug_t>::get_degree() {
    return neighbors.size();
}

template <typename aug_t>
bool ParallelUFOCluster<aug_t>::contracts() {
    if (!parent) return false;
    for (auto neighbor : neighbors)
        if (neighbor->parent == parent)
            return true;
    return false;
}

template <typename aug_t>
bool ParallelUFOCluster<aug_t>::atomic_contracts() {
    auto parent = AtomicLoad(&this->parent);
    if (!parent) return false;
    for (auto neighbor : neighbors)
        if (AtomicLoad(&neighbor->parent) == parent)
            return true;
    return false;
}

template <typename aug_t>
ParallelUFOCluster<aug_t>* ParallelUFOCluster<aug_t>::get_neighbor() {
    return *neighbors.begin();
}

template <typename aug_t>
ParallelUFOCluster<aug_t>* ParallelUFOCluster<aug_t>::get_other_neighbor(ParallelUFOCluster* c) {
    for (auto neighbor : neighbors)
        if (neighbor != c)
            return neighbor;
    return nullptr;
}

template <typename aug_t>
ParallelUFOCluster<aug_t>* ParallelUFOCluster<aug_t>::get_root() {
    ParallelUFOCluster* curr = this;
    while (curr->parent) curr = curr->parent;
    return curr;
}

}

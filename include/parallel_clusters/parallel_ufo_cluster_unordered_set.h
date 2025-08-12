#pragma once
#include "bridge.h"
#include "types.h"
#include "util.h"
#include <mutex>
#include <unordered_set>


namespace dgbs {


template <typename aug_t>
struct ParallelUFOClusterUSet {
  // Topology cluster data
  std::unordered_set<ParallelUFOClusterUSet*> neighbors;
  ParallelUFOClusterUSet* parent;
  ParallelUFOClusterUSet* partner;
  size_t granularity = 10000;
  // Constructor
  ParallelUFOClusterUSet() : parent(nullptr), partner(nullptr), neighbors() {};
  ParallelUFOClusterUSet(size_t _granularity) : parent(nullptr), partner(nullptr), neighbors() {granularity = _granularity;}
  // Functions
  std::mutex mtx;
  void insert_neighbor(ParallelUFOClusterUSet* c);
  void delete_neighbor(ParallelUFOClusterUSet* c);
  void insert_neighbors(parlay::sequence<ParallelUFOClusterUSet*>& cs);
  void delete_neighbors(parlay::sequence<ParallelUFOClusterUSet*>& cs);

  bool contains_neighbor(ParallelUFOClusterUSet* c);

  int get_degree();
  bool contracts();
  ParallelUFOClusterUSet* get_neighbor();
  ParallelUFOClusterUSet* get_other_neighbor(ParallelUFOClusterUSet* c);
  ParallelUFOClusterUSet* get_root();

  void print_neighbors();
};

// NOTE: Merge these insert and delete functions into 1 and use parfor with granularity.

template <typename aug_t>
void ParallelUFOClusterUSet<aug_t>::insert_neighbor(ParallelUFOClusterUSet* c) {
  mtx.lock();
  neighbors.insert(c);
  mtx.unlock();
}

template <typename aug_t>
void ParallelUFOClusterUSet<aug_t>::delete_neighbor(ParallelUFOClusterUSet* c) {
  mtx.lock();
  neighbors.erase(c);
  mtx.unlock();
}

template <typename aug_t>
void ParallelUFOClusterUSet<aug_t>::insert_neighbors(parlay::sequence<ParallelUFOClusterUSet*>& cs) {
  parlay::parallel_for(0, cs.size() - 1, [&] (int i){
    mtx.lock();
    neighbors.insert(cs[i]);
    mtx.unlock();
  },granularity);
}

template <typename aug_t>
void ParallelUFOClusterUSet<aug_t>::delete_neighbors(parlay::sequence<ParallelUFOClusterUSet*>& cs) {
  parlay::parallel_for(0, cs.size() - 1, [&] (int i){
    mtx.lock();
    neighbors.erase(cs[i]);
    mtx.unlock();
  }, granularity);

}


template <typename aug_t>
bool ParallelUFOClusterUSet<aug_t>::contains_neighbor(ParallelUFOClusterUSet* c) {
  return neighbors.contains(c);
}

template <typename aug_t>
int ParallelUFOClusterUSet<aug_t>::get_degree() {
  return neighbors.size();
}

template <typename aug_t>
bool ParallelUFOClusterUSet<aug_t>::contracts() {
  if (!parent) return false;
  for (auto neighbor : neighbors)
  if (neighbor->parent == parent)
    return true;
  return false;
}

template <typename aug_t>
ParallelUFOClusterUSet<aug_t>* ParallelUFOClusterUSet<aug_t>::get_neighbor() {
  return *neighbors.begin();
}

template <typename aug_t>
ParallelUFOClusterUSet<aug_t>* ParallelUFOClusterUSet<aug_t>::get_other_neighbor(ParallelUFOClusterUSet* c) {
  for (auto neighbor : neighbors)
  if (neighbor != c)
    return neighbor;
  return nullptr;
}

template <typename aug_t>
ParallelUFOClusterUSet<aug_t>* ParallelUFOClusterUSet<aug_t>::get_root() {
  ParallelUFOClusterUSet* curr = this;
  while (curr->parent) curr = curr->parent;
  return curr;
}

}

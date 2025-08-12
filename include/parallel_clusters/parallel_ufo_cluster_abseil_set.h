#pragma once
#include "absl/container/flat_hash_set.h"
#include "types.h"
#include "util.h"
#include <mutex>
#include <unordered_set>


namespace dgbs {

template <typename aug_t>
struct ParallelUFOClusterASet {
  // Topology cluster data
  absl::flat_hash_set<ParallelUFOClusterASet*> neighbors;
  ParallelUFOClusterASet* parent;
  ParallelUFOClusterASet* partner;
  size_t granularity = 10000;
  // Constructor
  ParallelUFOClusterASet() : parent(nullptr), partner(nullptr), neighbors() {};
  ParallelUFOClusterASet(size_t _granularity) : parent(nullptr), partner(nullptr), neighbors() {granularity = _granularity;}
  // Functions
  std::mutex mtx;
  void insert_neighbor(ParallelUFOClusterASet* c);
  void delete_neighbor(ParallelUFOClusterASet* c);
  void insert_neighbors(parlay::sequence<ParallelUFOClusterASet*>& cs);
  void delete_neighbors(parlay::sequence<ParallelUFOClusterASet*>& cs);

  bool contains_neighbor(ParallelUFOClusterASet* c);

  int get_degree();
  bool contracts();
  ParallelUFOClusterASet* get_neighbor();
  ParallelUFOClusterASet* get_other_neighbor(ParallelUFOClusterASet* c);
  ParallelUFOClusterASet* get_root();

  void print_neighbors();
};

// NOTE: Merge these insert and delete functions into 1 and use parfor with granularity.

template <typename aug_t>
void ParallelUFOClusterASet<aug_t>::insert_neighbor(ParallelUFOClusterASet* c) {
  mtx.lock();
  neighbors.insert(c);
  mtx.unlock();
}

template <typename aug_t>
void ParallelUFOClusterASet<aug_t>::delete_neighbor(ParallelUFOClusterASet* c) {
  mtx.lock();
  neighbors.erase(c);
  mtx.unlock();
}

template <typename aug_t>
void ParallelUFOClusterASet<aug_t>::insert_neighbors(parlay::sequence<ParallelUFOClusterASet*>& cs) {
  parlay::parallel_for(0, cs.size() - 1, [&] (int i){
    mtx.lock();
    neighbors.insert(cs[i]);
    mtx.unlock();
  }, granularity);

}

template <typename aug_t>
void ParallelUFOClusterASet<aug_t>::delete_neighbors(parlay::sequence<ParallelUFOClusterASet*>& cs) {
  parlay::parallel_for(0, cs.size() - 1, [&] (int i){
    mtx.lock();
    neighbors.erase(cs[i]);
    mtx.unlock();
  }, granularity);


}


template <typename aug_t>
bool ParallelUFOClusterASet<aug_t>::contains_neighbor(ParallelUFOClusterASet* c) {
  return neighbors.contains(c);
}

template <typename aug_t>
int ParallelUFOClusterASet<aug_t>::get_degree() {
  return neighbors.size();
}

template <typename aug_t>
bool ParallelUFOClusterASet<aug_t>::contracts() {
  if (!parent) return false;
  for (auto neighbor : neighbors)
  if (neighbor->parent == parent)
    return true;
  return false;
}

template <typename aug_t>
ParallelUFOClusterASet<aug_t>* ParallelUFOClusterASet<aug_t>::get_neighbor() {
  return *neighbors.begin();
}

template <typename aug_t>
ParallelUFOClusterASet<aug_t>* ParallelUFOClusterASet<aug_t>::get_other_neighbor(ParallelUFOClusterASet* c) {
  for (auto neighbor : neighbors)
  if (neighbor != c)
    return neighbor;
  return nullptr;
}

template <typename aug_t>
ParallelUFOClusterASet<aug_t>* ParallelUFOClusterASet<aug_t>::get_root() {
  ParallelUFOClusterASet* curr = this;
  while (curr->parent) curr = curr->parent;
  return curr;
}

}


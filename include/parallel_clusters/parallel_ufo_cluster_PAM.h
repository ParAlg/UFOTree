#pragma once
#include "absl/container/flat_hash_set.h"
#include "bridge.h"
#include "pam/utils.h"
#include "types.h"
#include "util.h"
#include <mutex>
#include <unordered_set>
#include "pam/pam.h"


namespace dgbs {

template<typename aug_t>
struct ParallelUFOClusterPAM;

struct entry{
  using key_t = ParallelUFOClusterPAM<int>*; 
  static inline bool comp(key_t a, key_t b) { return a < b;}

};

template <typename aug_t>
struct ParallelUFOClusterPAM {

  // Topology cluster data
  pam_set<entry> neighbors;
  ParallelUFOClusterPAM* parent;
  ParallelUFOClusterPAM* partner;
  size_t granularity = 1; // NOTE: Does not do anything, determined by PAM itself. Added for consistency between clusters.
  // Constructor
  ParallelUFOClusterPAM() : parent(nullptr), partner(nullptr), neighbors() {};
  ParallelUFOClusterPAM(size_t _granularity) : parent(nullptr), partner(nullptr), neighbors() {
    granularity = _granularity;
    //utils::node_limit = _granularity;
  }
  // Functions
  std::mutex mtx;
  void insert_neighbor(ParallelUFOClusterPAM* c);
  void delete_neighbor(ParallelUFOClusterPAM* c);
  void insert_neighbors(parlay::sequence<ParallelUFOClusterPAM*>& cs);
  void delete_neighbors(parlay::sequence<ParallelUFOClusterPAM*>& cs);

  bool contains_neighbor(ParallelUFOClusterPAM* c);

  int get_degree();
  bool contracts();
  ParallelUFOClusterPAM* get_neighbor();
  ParallelUFOClusterPAM* get_other_neighbor(ParallelUFOClusterPAM* c);
  ParallelUFOClusterPAM* get_root();

  void print_neighbors();
};

using pam_set = pam_set<entry>;
// NOTE: Merge these insert and delete functions into 1 and use parfor with granularity.
template <typename aug_t>
void ParallelUFOClusterPAM<aug_t>::insert_neighbor(ParallelUFOClusterPAM* c) {
  pam_set::insert(neighbors,c);
}

template <typename aug_t>
void ParallelUFOClusterPAM<aug_t>::delete_neighbor(ParallelUFOClusterPAM* c) {
  pam_set::remove(neighbors,c);
}

template <typename aug_t>
void ParallelUFOClusterPAM<aug_t>::insert_neighbors(parlay::sequence<ParallelUFOClusterPAM*>& cs) {
  pam_set::multi_insert(neighbors, cs);
}

template <typename aug_t>
void ParallelUFOClusterPAM<aug_t>::delete_neighbors(parlay::sequence<ParallelUFOClusterPAM*>& cs) {
    pam_set::multi_delete(neighbors, cs);
}


template <typename aug_t>
bool ParallelUFOClusterPAM<aug_t>::contains_neighbor(ParallelUFOClusterPAM* c) {
  return neighbors.contains(c);
}

template <typename aug_t>
int ParallelUFOClusterPAM<aug_t>::get_degree() {
  return neighbors.size();
}

template <typename aug_t>
bool ParallelUFOClusterPAM<aug_t>::contracts() {
  if (!parent) return false;
  /*for (auto neighbor : neighbors)
  if (neighbor->parent == parent)
    return true;i*/
  return false;
}

template <typename aug_t>
ParallelUFOClusterPAM<aug_t>* ParallelUFOClusterPAM<aug_t>::get_neighbor() {
  return nullptr; 
}

template <typename aug_t>
ParallelUFOClusterPAM<aug_t>* ParallelUFOClusterPAM<aug_t>::get_other_neighbor(ParallelUFOClusterPAM* c) {
  /*for (auto neighbor : neighbors)
  if (neighbor != c)
    return neighbor;*/
  return nullptr;
}

template <typename aug_t>
ParallelUFOClusterPAM<aug_t>* ParallelUFOClusterPAM<aug_t>::get_root() {
  ParallelUFOClusterPAM* curr = this;
  while (curr->parent) curr = curr->parent;
  return curr;
}

}


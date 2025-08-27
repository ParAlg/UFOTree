#pragma once
#include <parlay/alloc.h>
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include "bridge.h"
#include "sparse_table.h"
#include "types.h"
#include "util.h"

using namespace parlay;

namespace dgbs {
template <typename aug_t>
struct ParallelTopologyCluster {
  // Topology cluster data
  ParallelTopologyCluster<aug_t>* neighbors[3];
  ParallelTopologyCluster<aug_t>* parent;
  ParallelTopologyCluster<aug_t>* partner;

  aug_t edge_values[3];   // Only for path queries
  aug_t value;            // Stores subtree values or cluster path values
  // Some fields for asynchronous operations
  uint32_t priority;
  bool del;

  // Constructor
  ParallelTopologyCluster(){};
  ParallelTopologyCluster(uint32_t salt, aug_t value = 1)
      : neighbors(),
        edge_values(),
        value(value),
        parent(nullptr),
        del(false),
        partner(nullptr) {
    priority = hash32(salt);
  };
  // Helper functions
  int get_degree();
  bool contracts();
  bool contains_neighbor(ParallelTopologyCluster<aug_t>* c);
  void insert_neighbor(ParallelTopologyCluster<aug_t>* c, aug_t value);
  void remove_neighbor(ParallelTopologyCluster<aug_t>* c);
  ParallelTopologyCluster<aug_t>* get_neighbor(int index = 0);
  ParallelTopologyCluster<aug_t>* get_other_neighbor(
     ParallelTopologyCluster<aug_t>* first_neighbor);
  ParallelTopologyCluster<aug_t>* get_root();

  // Try to delete this cluster non-atomically.
  bool try_del();
  // Try to delete this cluster *atomically* (i.e., by CAS'ing the del
  // field).
  bool try_del_atomic();
};

template <typename aug_t>
int ParallelTopologyCluster<aug_t>::get_degree() {
  int deg = 0;
  for (auto neighbor : this->neighbors)
    if (neighbor)
      deg++;
  return deg;
}

// Helper function which returns whether this cluster combines with another
// cluster.
template <typename aug_t>
bool ParallelTopologyCluster<aug_t>::contracts() {
  if (!parent)
    return false;
  bool contracts = false;
  for (auto neighbor : this->neighbors)
    if (neighbor && neighbor->parent == this->parent)
      contracts = true;
  return contracts;
}

template <typename aug_t>
bool ParallelTopologyCluster<aug_t>::contains_neighbor(
   ParallelTopologyCluster<aug_t>* c) {
  for (int i = 0; i < 3; ++i)
    if (this->neighbors[i] == c)
      return true;
  return false;
}

template <typename aug_t>
void ParallelTopologyCluster<aug_t>::insert_neighbor(
   ParallelTopologyCluster<aug_t>* c, aug_t value) {
  if (this->contains_neighbor(c))
    return;
  for (int i = 0; i < 3; ++i) {
    if (gbbs::CAS(&this->neighbors[i], (ParallelTopologyCluster<aug_t>*)nullptr, c)) {
      this->edge_values[i] = value;
      return;
    } else if (this->neighbors[i] == c)
      return;
  }
  std::cerr << "No space to insert neighbor." << std::endl;
  std::abort();
}

template <typename aug_t>
void ParallelTopologyCluster<aug_t>::remove_neighbor(
   ParallelTopologyCluster<aug_t>* c) {
  for (int i = 0; i < 3; ++i) {
    if (this->neighbors[i] == c) {
      this->neighbors[i] = nullptr;
      return;
    }
  }
  std::cerr << "Neighbor to delete not found." << std::endl;
  std::abort();
}

template <typename aug_t>
ParallelTopologyCluster<aug_t>* ParallelTopologyCluster<aug_t>::get_neighbor(
   int index) {
  assert(get_degree() > index);
  int neighbors_seen = 0;
  for (auto neighbor : neighbors) {
    if (neighbor) {
      if (neighbors_seen == index)
        return neighbor;
      neighbors_seen++;
    }
  }
  return nullptr;
}

template <typename aug_t>
ParallelTopologyCluster<aug_t>*
ParallelTopologyCluster<aug_t>::get_other_neighbor(
   ParallelTopologyCluster<aug_t>* first_neighbor) {
  assert(contains_neighbor(first_neighbor));
  for (auto neighbor : neighbors)
    if (neighbor && neighbor != first_neighbor)
      return neighbor;
  return nullptr;
}

template <typename aug_t>
ParallelTopologyCluster<aug_t>* ParallelTopologyCluster<aug_t>::get_root() {
  ParallelTopologyCluster<aug_t>* curr = this;
  while (curr->parent)
    curr = curr->parent;
  return curr;
}

template <typename aug_t>
bool ParallelTopologyCluster<aug_t>::try_del() {
  if (!del) {
    del = true;
    return true;
  }
  return false;
}

template <typename aug_t>
bool ParallelTopologyCluster<aug_t>::try_del_atomic() {
  if (!del) {
    bool ret = gbbs::CAS(&del, false, true);
    return ret;
  }
  return false;
}
}
#pragma once

// This is a naive augmented skip list where a batch of k joins/splits just does each
// sequentially, updating the augmented value after each join/split. Work bound
// is O(k log n) for k operations over n elements.
//
// The augmentation is hardcoded to the sum function with the value 1 assigned
// to each element. As such, `GetSum()` returns the size of the list.

#include <algorithm>
#include <cassert>
#include <tuple>

#include <parett/sequence/skip_list/include/skip_list_base.hpp>

using std::pair;
using std::tie;


namespace skip_list {

class AugmentedElement : private ElementBase<AugmentedElement> {
  friend class ElementBase<AugmentedElement>;
 public:
  // See comments on [ElementBase<>]
  AugmentedElement();
  AugmentedElement(size_t random_int);

  ~AugmentedElement();

  // Get result of applying augmentation function over the whole list.
  // May run concurrently with other [GetAggregate] calls.
  int GetSum();

  // For each pair [<left, right>] in [joins], [left] must be the last node in
  // its list, and [right] must be the first node of in its list. Each [left]
  // must be unique, and each [right] must be unique.
  static void BatchJoin(
      std::pair<AugmentedElement*, AugmentedElement*>* joins, int len);

  static void BatchSplit(AugmentedElement** splits, int len);

  using ElementBase<AugmentedElement>::FindRepresentative;
  using ElementBase<AugmentedElement>::GetNextElement;
  using ElementBase<AugmentedElement>::GetPreviousElement;
  
  size_t calculate_size(){
    size_t memory = ElementBase<AugmentedElement>::calculate_size();
    int* curr = vals;
    while(curr){
      memory += sizeof(curr);
      curr++;
    }
    return memory;
  }
 private:
  int* vals;

  void AllocateVals();

  // return left parent & sum from left parent's direct child (inclusive) to
  // input node (inclusive)
  std::pair<AugmentedElement*, int> FindLeftParentAndSum(int level);
  // return right parent & sum from input node (inclusive) to right parent's //
  // direct child (exclusive)
  std::pair<AugmentedElement*, int> FindRightParentAndSum(int level);
};


void AugmentedElement::AllocateVals() {
  vals = new int[height];
  for (int i = 0; i < height; i++) {
    vals[i] = 1;
  }
}

AugmentedElement::AugmentedElement() :
  ElementBase<AugmentedElement>() {
  AllocateVals();
}

AugmentedElement::AugmentedElement(size_t random_int) :
  ElementBase<AugmentedElement>(random_int) {
  AllocateVals();
}

AugmentedElement::~AugmentedElement() {
  delete[] vals;
}

int AugmentedElement::GetSum() {
  AugmentedElement* root = FindRepresentative();
  int level = root->height - 1;
  int sum = root->vals[level];
  AugmentedElement* curr = root->neighbors[level].next;
  while (curr != nullptr && curr != root) {
    sum += curr->vals[level];
    curr = curr->neighbors[level].next;
  }
  if (curr == nullptr) {
    // list is not circular; need to traverse backwards to beginning of list
    curr = root;
    while (true) {
      while (level >= 0 && curr->neighbors[level].prev == nullptr) {
        level--;
      }
      if (level < 0) {
        break;
      }
      while (curr->neighbors[level].prev != nullptr) {
        curr = curr->neighbors[level].prev;
        sum += curr->vals[level];
      }
    }
  }
  return sum;
}

pair<AugmentedElement*, int> AugmentedElement::FindLeftParentAndSum(
    int level) {
  int sum = 0;
  AugmentedElement* current_node = this;
  AugmentedElement* start_node = current_node;
  do {
    sum += current_node->vals[level];
    if (current_node->height > level + 1) {
      return make_pair(current_node, sum);
    }
    current_node = current_node->neighbors[level].prev;
  } while (current_node != nullptr && current_node != start_node);
  return make_pair(nullptr, sum);
}

pair<AugmentedElement*, int> AugmentedElement::FindRightParentAndSum(
    int level) {
  int sum = 0;
  AugmentedElement* current_node = this;
  AugmentedElement* start_node = current_node;
  do {
    if (current_node->height > level + 1) {
      return make_pair(current_node, sum);
    }
    sum += current_node->vals[level];
    current_node = current_node->neighbors[level].next;
  } while (current_node != nullptr && current_node != start_node);
  return make_pair(nullptr, sum);
}

void AugmentedElement::BatchJoin(
    pair<AugmentedElement*, AugmentedElement*>* joins, int len) {
  for (int i = 0; i < len; i++) {
    AugmentedElement* left, * right;
    tie(left, right) = joins[i];

    int level = 0;
    int sum = left->vals[0];
    while (left != nullptr && right != nullptr) {
      left->neighbors[level].next = right;
      right->neighbors[level].prev = left;
      left->vals[level] = sum;
      int sum_left, sum_right;
      tie(left, sum_left) = left->FindLeftParentAndSum(level);
      tie(right, sum_right) = right->FindRightParentAndSum(level);
      sum = sum_left + sum_right;
      level++;
    }
    while (left != nullptr) {
      left->vals[level] = sum;
      tie(left, sum) = left->FindLeftParentAndSum(level);
      level++;
    }
  }
}

void AugmentedElement::BatchSplit(AugmentedElement** splits, int len) {
  for (int i = 0; i < len; i++) {
    AugmentedElement* current = splits[i];
    int level = 0;
    int sum = current->vals[0];
    AugmentedElement* next;
    while (current != nullptr &&
        (next = current->neighbors[level].next) != nullptr) {
      current->neighbors[level].next = nullptr;
      next->neighbors[level].prev = nullptr;
      current->vals[level] = sum;
      tie(current, sum) = current->FindLeftParentAndSum(level);
      level++;
    }
    while (current != nullptr) {
      current->vals[level] = sum;
      tie(current, sum) = current->FindLeftParentAndSum(level);
      level++;
    }
  }
}

} // namespace skip_list

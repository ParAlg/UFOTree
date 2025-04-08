#pragma once

#include <utility>
#include <tuple>
#include <parett/utilities/blockRadixSort.h>
#include <parett/utilities/random.h>
#include "types.h"

using std::pair;
using namespace dgbs;


namespace treap {

template<typename T>
class Node {
 public:
  Node();
  explicit Node(unsigned random_int);

  Node* GetRoot() const;

  // SplitRights right after this node
  std::pair<Node*, Node*> SplitRight();
  std::pair<Node*, Node*> SplitLeft();
  std::pair<Node*, Node*> Split() { return SplitRight(); };
  std::pair<Node*, Node*> SplitAround();
  // Join tree containing lesser to tree containing greater and return root of
  // resulting tree.
  static Node* Join(Node* lesser, Node* greater);

  static std::function<T(T, T)> aggregate_function;

 private:
 static Node* JoinRoots(Node* lesser, Node* greater);
  void AssignChild(int i, Node* v);
  void RemoveChild(int i);
  void RecomputeAggregate();

  [[no_unique_address]] T value;
  [[no_unique_address]] T aggregate;
  Node* parent_;
  Node* child_[2];
  uint32_t priority_;
};


namespace {
  pbbs::random default_randomness;
}  // namespace

template<typename T>
Node<T>::Node(unsigned random_int)
  : parent_(nullptr)
  , child_{nullptr, nullptr}
  , priority_(random_int) {}

template<typename T>
Node<T>::Node() : Node(default_randomness.rand()) {
  default_randomness = default_randomness.next();
}

template<typename T>
std::function<T(T,T)> Node<T>::aggregate_function = [] (T x, T y) { return x; };

template<typename T>
inline void Node<T>::AssignChild(int i, Node* v) {
  if (v != nullptr)
    v->parent_ = this;
  child_[i] = v;
  if constexpr (!std::is_same<T, empty_t>::value)
    RecomputeAggregate();
}

template<typename T>
inline void Node<T>::RemoveChild(int i) {
  if (child_[i])
    child_[i]->parent_ = nullptr;
  child_[i] = nullptr;
  if constexpr (!std::is_same<T, empty_t>::value)
    RecomputeAggregate();
}

template<typename T>
inline void Node<T>::RecomputeAggregate() {
  aggregate = value;
  if (child_[0])
    aggregate = aggregate_function(aggregate, child_[0]->aggregate);
  if (child_[1])
    aggregate = aggregate_function(aggregate, child_[1]->aggregate);
}

template<typename T>
Node<T>* Node<T>::GetRoot() const {
  const Node* current = this;
  while (current->parent_ != nullptr) {
    current = current->parent_;
  }
  return const_cast<Node*>(current);
}

template<typename T>
pair<Node<T>*, Node<T>*> Node<T>::SplitRight() {
  Node* lesser = nullptr;
  Node* greater = child_[1];
  RemoveChild(1);

  Node* current = this;
  bool traversed_up_from_right = 1;
  bool next_direction;
  while (current != nullptr) {
    Node* p = current->parent_;
    if (p != nullptr) {
      next_direction = p->child_[1] == current;
      p->RemoveChild(next_direction);
    }
    if (traversed_up_from_right) {
      lesser = Join(current, lesser);
    } else {
      greater = Join(greater, current);
    }

    traversed_up_from_right = next_direction;
    current = p;
  }
  return {lesser, greater};
}

template<typename T>
pair<Node<T>*, Node<T>*> Node<T>::SplitLeft() {
  Node* lesser = child_[0];
  Node* greater = nullptr;
  RemoveChild(0);

  Node* current = this;
  bool traversed_up_from_right = 1;
  bool next_direction;
  while (current != nullptr) {
    Node* p = current->parent_;
    if (p != nullptr) {
      next_direction = p->child_[1] == current;
      p->RemoveChild(next_direction);
    }
    if (traversed_up_from_right) {
      lesser = Join(current, lesser);
    } else {
      greater = Join(greater, current);
    }

    traversed_up_from_right = next_direction;
    current = p;
  }
  return {lesser, greater};
}

template<typename T>
pair<Node<T>*, Node<T>*> Node<T>::SplitAround() {
  Node* lesser = child_[0];
  Node* greater = child_[1];
  RemoveChild(0);
  RemoveChild(1);

  Node* current = this;
  bool traversed_up_from_right;
  bool next_direction;
  while (current != nullptr) {
    Node* p = current->parent_;
    if (p != nullptr) {
      next_direction = p->child_[1] == current;
      p->RemoveChild(next_direction);
    }
    if (current != this && traversed_up_from_right) {
      lesser = Join(current, lesser);
    } else {
      greater = Join(greater, current);
    }

    traversed_up_from_right = next_direction;
    current = p;
  }
  return {lesser, greater};
}

template<typename T>
Node<T>* Node<T>::JoinRoots(Node* lesser, Node* greater) {
  if (lesser == nullptr) {
    return greater;
  } else if (greater == nullptr) {
    return lesser;
  }

  if (lesser->priority_ > greater->priority_) {
    lesser->AssignChild(1, JoinRoots(lesser->child_[1], greater));
    return lesser;
  } else {
    greater->AssignChild(0, JoinRoots(lesser, greater->child_[0]));
    return greater;
  }
}

template<typename T>
Node<T>* Node<T>::Join(Node* lesser, Node* greater) {
  return JoinRoots(
      lesser == nullptr ? nullptr : lesser->GetRoot(),
      greater == nullptr ? nullptr : greater->GetRoot());
}

}  // namespace treap

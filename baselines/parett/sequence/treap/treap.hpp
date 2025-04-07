#pragma once

#include <utility>
#include <tuple>
#include <parett/utilities/blockRadixSort.h>
#include <parett/utilities/random.h>


using std::pair;

namespace treap {

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

 private:
  void AssignChild(int i, Node* v);
  static Node* JoinRoots(Node* lesser, Node* greater);

  Node* parent_;
  Node* child_[2];
  unsigned priority_;
};


namespace {
  pbbs::random default_randomness;
}  // namespace

Node::Node(unsigned random_int)
  : parent_(nullptr)
  , child_{nullptr, nullptr}
  , priority_(random_int) {}

Node::Node() : Node(default_randomness.rand()) {
  default_randomness = default_randomness.next();
}

void Node::AssignChild(int i, Node* v) {
  if (v != nullptr) {
    v->parent_ = this;
  }
  child_[i] = v;
}

Node* Node::GetRoot() const {
  const Node* current = this;
  while (current->parent_ != nullptr) {
    current = current->parent_;
  }
  return const_cast<Node*>(current);
}

pair<Node*, Node*> Node::SplitRight() {
  Node* lesser = nullptr;
  Node* greater = child_[1];
  if (child_[1] != nullptr) {
    child_[1]->parent_ = nullptr;
    AssignChild(1, nullptr);
  }

  Node* current = this;
  bool traversed_up_from_right = 1;
  bool next_direction;
  while (current != nullptr) {
    Node* p = current->parent_;
    if (p != nullptr) {
      next_direction = p->child_[1] == current;
      p->AssignChild(next_direction, nullptr);
      current->parent_ = nullptr;
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

pair<Node*, Node*> Node::SplitLeft() {
  Node* lesser = child_[0];
  Node* greater = nullptr;
  if (child_[0] != nullptr) {
    child_[0]->parent_ = nullptr;
    AssignChild(0, nullptr);
  }

  Node* current = this;
  bool traversed_up_from_right = 1;
  bool next_direction;
  while (current != nullptr) {
    Node* p = current->parent_;
    if (p != nullptr) {
      next_direction = p->child_[1] == current;
      p->AssignChild(next_direction, nullptr);
      current->parent_ = nullptr;
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

pair<Node*, Node*> Node::SplitAround() {
  Node* lesser = child_[0];
  Node* greater = child_[1];
  if (child_[0] != nullptr) {
    child_[0]->parent_ = nullptr;
    AssignChild(0, nullptr);
  }
  if (child_[1] != nullptr) {
    child_[1]->parent_ = nullptr;
    AssignChild(1, nullptr);
  }

  Node* current = this;
  bool traversed_up_from_right;
  bool next_direction;
  while (current != nullptr) {
    Node* p = current->parent_;
    if (p != nullptr) {
      next_direction = p->child_[1] == current;
      p->AssignChild(next_direction, nullptr);
      current->parent_ = nullptr;
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

Node* Node::JoinRoots(Node* lesser, Node* greater) {
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

Node* Node::Join(Node* lesser, Node* greater) {
  return JoinRoots(
      lesser == nullptr ? nullptr : lesser->GetRoot(),
      greater == nullptr ? nullptr : greater->GetRoot());
}

}  // namespace treap

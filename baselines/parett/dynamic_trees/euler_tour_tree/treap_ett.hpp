#pragma once

#include <algorithm>
#include <unordered_map>

#include <parett/sequence/treap/treap.hpp>
#include <parett/utilities/hash_pair.hpp>


namespace dgbs {

namespace treap {

using std::pair;

template<typename T>
class EulerTourTree {
using Node = treap::Node<T>;

public:
  EulerTourTree(vertex_t _num_verts);
  ~EulerTourTree();
  EulerTourTree(const EulerTourTree& other);
  EulerTourTree& operator=(const EulerTourTree& other);

  bool IsConnected(vertex_t u, vertex_t v);
  void Link(vertex_t u, vertex_t v);
  void LinkInner(vertex_t u, vertex_t v);
  void LinkOuter(vertex_t u, vertex_t v);
  void Cut(vertex_t u, vertex_t v);

  bool connected(vertex_t u, vertex_t v) { return IsConnected(u,v); }
  void link(vertex_t u, vertex_t v) { Link(u,v); }
  void cut(vertex_t u, vertex_t v) { Cut(u,v); }

  T GetValue(vertex_t v);
  T GetNodeAggregate(Node* node);
  T GetComponentAggregate(vertex_t v);
  T GetSubtreeAggregate(vertex_t v, vertex_t p);
  
  void UpdateValue(vertex_t v, T value);
  void UpdateWithFunction(vertex_t v, std::function<bool(T&)> f);

  Node* GetVertexNode(vertex_t v);
  Node* GetEdgeNode(vertex_t u, vertex_t v);
  vertex_t NodeToVertex(Node* node);
  std::unordered_map<std::pair<vertex_t, vertex_t>, Node*, HashIntPairStruct>& GetEdgeMap() { return edges; }

private:
  vertex_t num_verts;
  Node* verts;
  std::unordered_map<std::pair<vertex_t, vertex_t>, Node*, HashIntPairStruct> edges;
  std::vector<Node*> node_pool;
};


template<typename T>
EulerTourTree<T>::EulerTourTree(vertex_t _num_verts) : num_verts(_num_verts) {
  verts = new Node[num_verts];
  edges.reserve(2 * (num_verts - 1));
  for (vertex_t i = 0; i < 2 * (num_verts - 1); i++)
    node_pool.push_back(new Node());
}

template<typename T>
EulerTourTree<T>::~EulerTourTree() {
  delete[] verts;
  for (auto it : edges)
    delete it.second;
  for (Node* node : node_pool)
    delete node;
}

template<typename T>
bool EulerTourTree<T>::IsConnected(vertex_t u, vertex_t v) {
  return verts[u].GetRoot() == verts[v].GetRoot();
}

template<typename T>
void EulerTourTree<T>::Link(vertex_t u, vertex_t v) {
  Node* u_left, * u_right, * v_left, * v_right;
  Node* uv = node_pool.back();
  node_pool.pop_back();
  Node* vu = node_pool.back();
  node_pool.pop_back();
  edges[std::make_pair(u, v)] = uv;
  edges[std::make_pair(v, u)] = vu;
  std::tie(u_left, u_right) = verts[u].SplitRight();
  std::tie(v_left, v_right) = verts[v].SplitRight();
  Node* root_left = Node::Join(u_left, uv);
  root_left = Node::Join(root_left, v_right);
  Node* root_right = Node::Join(v_left, vu);
  root_right = Node::Join(root_right, u_right);
  Node::Join(root_left, root_right);
}

template<typename T>
void EulerTourTree<T>::LinkInner(vertex_t u, vertex_t v) {
  Node* u_left, * u_right, * v_left, * v_right;
  Node* uv = node_pool.back();
  node_pool.pop_back();
  Node* vu = node_pool.back();
  node_pool.pop_back();
  edges[std::make_pair(u, v)] = uv;
  edges[std::make_pair(v, u)] = vu;
  std::tie(u_left, u_right) = verts[u].SplitRight();
  std::tie(v_left, v_right) = verts[v].SplitLeft();
  Node* root_left = Node::Join(u_left, uv);
  root_left = Node::Join(root_left, v_right);
  Node* root_right = Node::Join(v_left, vu);
  root_right = Node::Join(root_right, u_right);
  Node::Join(root_left, root_right);
}

template<typename T>
void EulerTourTree<T>::LinkOuter(vertex_t u, vertex_t v) {
  Node* u_left, * u_right, * v_left, * v_right;
  Node* uv = node_pool.back();
  node_pool.pop_back();
  Node* vu = node_pool.back();
  node_pool.pop_back();
  edges[std::make_pair(u, v)] = uv;
  edges[std::make_pair(v, u)] = vu;
  std::tie(u_left, u_right) = verts[u].SplitLeft();
  std::tie(v_left, v_right) = verts[v].SplitLeft();
  Node* root_left = Node::Join(u_left, uv);
  root_left = Node::Join(root_left, v_right);
  Node* root_right = Node::Join(v_left, vu);
  root_right = Node::Join(root_right, u_right);
  Node::Join(root_left, root_right);
}

template<typename T>
void EulerTourTree<T>::Cut(vertex_t u, vertex_t v) {
  Node* uv_left, * uv_right, * vu_left, * vu_right;
  auto uv_it = edges.find(std::make_pair(u, v));
  auto vu_it = edges.find(std::make_pair(v, u));
  std::tie(uv_left, uv_right) = uv_it->second->SplitAround();
  bool uv_vu_in_order = !uv_left || (uv_left->GetRoot() != vu_it->second->GetRoot());
  std::tie(vu_left, vu_right) = vu_it->second->SplitAround();
  node_pool.push_back(uv_it->second);
  node_pool.push_back(vu_it->second);
  edges.erase(uv_it);
  edges.erase(vu_it);
  if (uv_vu_in_order) {
    Node::Join(uv_left, vu_right);
  } else {
    Node::Join(vu_left, uv_right);
  }
}

template<typename T>
T EulerTourTree<T>::GetValue(vertex_t v) {
  return verts[v].value;
}

template<typename T>
T EulerTourTree<T>::GetNodeAggregate(Node* node) {
  return node->aggregate;
}

template<typename T>
T EulerTourTree<T>::GetComponentAggregate(vertex_t v) {
  return verts[v].GetRoot()->aggregate;
}

template<typename T>
T EulerTourTree<T>::GetSubtreeAggregate(vertex_t v, vertex_t p) {
  return Node::empty_augment;
}


template<typename T>
void EulerTourTree<T>::UpdateValue(vertex_t v, T value) {
  verts[v].value = value;
  Node* curr = &verts[v];
  while (curr) {
    curr->RecomputeAggregate();
    curr = curr->parent_;
  }
}

template<typename T>
void EulerTourTree<T>::UpdateWithFunction(vertex_t v, std::function<bool(T&)> f) {
  f(verts[v].value);
  Node* curr = &verts[v];
  while (curr && f(curr->aggregate)) {
    curr = curr->parent_;
  }
}

template<typename T>
treap::Node<T>* EulerTourTree<T>::GetVertexNode(vertex_t v) {
  return &verts[v];
}

template<typename T>
treap::Node<T>* EulerTourTree<T>::GetEdgeNode(vertex_t u, vertex_t v) {
  auto it = edges.find(std::make_pair(u,v));
  if (it == edges.end()) return nullptr;
  return it->second;
}

template<typename T>
vertex_t EulerTourTree<T>::NodeToVertex(Node* node) {
  return node - &verts[0];
}

} //namespace treap

}

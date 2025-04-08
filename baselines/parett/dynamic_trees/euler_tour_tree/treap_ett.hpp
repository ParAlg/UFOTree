#pragma once

#include <algorithm>
#include <unordered_map>

#include <parett/sequence/treap/treap.hpp>
#include <parett/utilities/hash_pair.hpp>


namespace dgbs {

namespace treap_ett {

using std::pair;

template<typename T>
class EulerTourTree {
using Node = treap::Node<T>;

public:
  EulerTourTree(int _num_verts);
  ~EulerTourTree();
  EulerTourTree(const EulerTourTree& other);
  EulerTourTree& operator=(const EulerTourTree& other);

  bool IsConnected(int u, int v);
  void Link(int u, int v);
  void Cut(int u, int v);

  bool connected(int u, int v){return IsConnected(u,v);}
  void link(int u, int v){Link(u,v);}
  void cut(int u, int v){ Cut(u,v);}

private:
  int num_verts;
  Node* verts;
  std::unordered_map<std::pair<int, int>, Node*, HashIntPairStruct> edges;
  std::vector<Node*> node_pool;
};


template<typename T>
EulerTourTree<T>::EulerTourTree(int _num_verts) : num_verts(_num_verts) {
  verts = new Node[num_verts];
  edges.reserve(2 * (num_verts - 1));
  for (int i = 0; i < 2 * (num_verts - 1); i++)
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
bool EulerTourTree<T>::IsConnected(int u, int v) {
  return verts[u].GetRoot() == verts[v].GetRoot();
}

template<typename T>
void EulerTourTree<T>::Link(int u, int v) {
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
void EulerTourTree<T>::Cut(int u, int v) {
  Node* uv_left, * uv_right, * vu_left, * vu_right;
  auto uv_it = edges.find(std::make_pair(u, v));
  auto vu_it = edges.find(std::make_pair(v, u));
  std::tie(uv_left, uv_right) = uv_it->second->SplitAround();
  std::tie(vu_left, vu_right) = vu_it->second->SplitAround();
  node_pool.push_back(uv_it->second);
  node_pool.push_back(vu_it->second);
  edges.erase(uv_it);
  edges.erase(vu_it);
  if (uv_left)  uv_left  = uv_left->GetRoot();
  if (uv_right) uv_right = uv_right->GetRoot();
  if (vu_left)  vu_left  = vu_left->GetRoot();
  if (vu_right) vu_right = vu_right->GetRoot();
  if (uv_left != vu_right) {
    Node::Join(uv_left, vu_right);
  } else {
    Node::Join(vu_left, uv_right);
  }
}

} //namespace treap_ett

}

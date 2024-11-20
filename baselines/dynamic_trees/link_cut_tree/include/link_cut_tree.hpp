#pragma once

#include <algorithm>
#include "types.h"

namespace link_cut_tree {

class Node;

class LinkCutTree {
 public:
  LinkCutTree(int _num_verts);
  ~LinkCutTree();
  
  void link(vertex_t u, vertex_t v);
  void cut(vertex_t u, vertex_t v);
  bool connected(vertex_t u, vertex_t v);

  bool* BatchConnected(std::pair<int, int>* queries, int len);
  // Inserting all links in [links] must keep the graph acylic.
  void BatchLink(std::pair<int, int>* links, int len);
  // All edges in [cuts] must be in the graph, and no edges may be repeated.
  void BatchCut(std::pair<int, int>* cuts, int len);

 private:
  Node* verts;
  int num_verts;
};


// ============= Link Cut Tree With Integer Path Queries Below This Point ===============

class NodeInt;

class LinkCutTreeInt {
 public:
  LinkCutTreeInt(int _num_verts);
  ~LinkCutTreeInt();
  
  void link(vertex_t u, vertex_t v, int weight = 0);
  void cut(vertex_t u, vertex_t v);
  bool connected(vertex_t u, vertex_t v);

  bool* BatchConnected(std::pair<int, int>* queries, int len);
  // Inserting all links in [links] must keep the graph acylic.
  void BatchLink(std::pair<int, int>* links, int len);
  // All edges in [cuts] must be in the graph, and no edges may be repeated.
  void BatchCut(std::pair<int, int>* cuts, int len);

 private:
  NodeInt* verts;
  int num_verts;
};

} // namespace link_cut_tree

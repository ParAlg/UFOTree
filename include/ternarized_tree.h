#pragma once
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <unordered_map>
#include <queue>
#include <cassert>
#include "util.h"

struct pair_hash{ 
  template<typename U, typename V>
  std::size_t operator()(const std::pair<U,V> x) const{
    return std::hash<long long>()(((long long)x.first)^(((long long)x.second)<<32));
  }
};

template<typename DynamicTree, typename aug_t>
class TernarizedTree {
public:
  // Ternarized tree interface
  TernarizedTree(vertex_t n, QueryType q = PATH, 
                 std::function<aug_t(aug_t, aug_t)> f = [](aug_t x, aug_t y)->aug_t{return x + y;}, 
                 aug_t id = 0, aug_t dval = 0);

  void link(vertex_t u, vertex_t v, aug_t value = 1);
  void cut(vertex_t u, vertex_t v);
  bool connected(vertex_t u, vertex_t v);
  //private: - Making public for testing
  vertex_t determine_link_v(vertex_t u);
  vertex_t ternarize_vertex(vertex_t u);
  void delete_ternarized_vertex(vertex_t v);

  // Test helper functions
  bool is_valid_ternarized_tree();
  int get_length_of_chain(vertex_t v);
  bool vertex_on_chain(vertex_t start, vertex_t to_find);
  bool test_vertex_deleted(vertex_t u, vertex_t v);
  // Underlying dynamic tree data structure
  DynamicTree tree;
  // Ternarization book-keeping
  vertex_t n;
  aug_t id;
  std::unordered_map<std::pair<vertex_t, vertex_t>, std::pair<vertex_t, vertex_t>, pair_hash> edge_map;
  std::unordered_map<vertex_t, std::pair<vertex_t, vertex_t>> chain_map;
  std::vector<vertex_t> head_map;
  std::queue<vertex_t> free_ids;
};

template<typename DynamicTree, typename aug_t>
TernarizedTree<DynamicTree, aug_t>::TernarizedTree(vertex_t n, QueryType q, 
                                                                std::function<aug_t(aug_t, aug_t)> f, 
                                                                aug_t id, aug_t d) : 
  tree(2*n, q, f, id, d), n(n) {
  head_map.resize(n, MAX_VERTEX_T);
  for(int i = n; i < (2*n); i++) free_ids.push(i);
}


/*** HELPER METHODS ***/

template<typename DynamicTree, typename aug_t>
vertex_t TernarizedTree<DynamicTree, aug_t>::determine_link_v(vertex_t v){
  vertex_t prev_ternary_neighbor = MAX_VERTEX_T;
  if(chain_map.find(v) != chain_map.end()) prev_ternary_neighbor = chain_map[v].second;
  if(prev_ternary_neighbor == MAX_VERTEX_T){
    prev_ternary_neighbor = ternarize_vertex(v); 
  } else {
    tree.cut(v, prev_ternary_neighbor);
  }
  vertex_t new_ternarized_n = free_ids.front(); free_ids.pop();
  head_map[new_ternarized_n - n] = v;
  tree.link(v, new_ternarized_n, id);
  tree.link(prev_ternary_neighbor, new_ternarized_n, id);
  // Adjustments to the chain.
  chain_map[v].second = new_ternarized_n; 
  chain_map[new_ternarized_n].first = v; chain_map[new_ternarized_n].second = prev_ternary_neighbor;
  chain_map[prev_ternary_neighbor].first = new_ternarized_n;
  // This is vertex to link with however the other vertex is handled.
  return new_ternarized_n;
}

/* Given a degree 3 vertex, ternarize the vertex by replacing one of its edges to a
 * vertex with an edge to a ternarized vertex that as an edge to the original vertex
 * 
 * Ex: (d is the new ternarized vertex)
 *        a-u-b
 *          |
 *          c
 * becomes 
 *        a-u-b
 *          |
 *          d-c*/
template<typename DynamicTree, typename aug_t>
vertex_t TernarizedTree<DynamicTree, aug_t>::ternarize_vertex(vertex_t v){
  assert(tree.get_degree(v) == 3);
  auto to_cut_pair = tree.retrieve_v_to_del(v); // Retrieve the vertex, edge pair to add as the start of the ternarization chain.
  tree.cut(v, to_cut_pair.first); // cut edge between this vertex to add a dummy ternarization vertex.
  vertex_t new_ternarized_n = free_ids.front(); free_ids.pop();

  auto head = to_cut_pair.first >=n ? head_map[to_cut_pair.first - n] : to_cut_pair.first;
  if(head > v){ 
    edge_map[std::pair(v, head)].first = new_ternarized_n;
  } else {
    edge_map[std::pair(head, v)].second = new_ternarized_n;
  }
  tree.link(new_ternarized_n, to_cut_pair.first, to_cut_pair.second);
  head_map[new_ternarized_n - n] = v;
  chain_map[v] = std::pair(MAX_VERTEX_T, new_ternarized_n); // Head of new ternarized chain.
  chain_map[new_ternarized_n] = std::pair(v, MAX_VERTEX_T); // Add as tail of ternarized chain. 
  return new_ternarized_n;
}

template<typename DynamicTree, typename aug_t>
void TernarizedTree<DynamicTree, aug_t>::delete_ternarized_vertex(vertex_t v){
  // If vertex is degree 1, it is at the tail of the ternarization chain.
  if(tree.get_degree(v) == 1){
    vertex_t neighbor = chain_map[v].first; // Get prev reference of vertex.
    tree.cut(v, neighbor); // Delete this edge.
    // Adjustments to chain.
    chain_map.erase(v);
    chain_map[neighbor].second = MAX_VERTEX_T;
    if(neighbor == head_map[v-n]) chain_map.erase(neighbor);
  } else {
    vertex_t neighbor1 = chain_map[v].first, neighbor2 = chain_map[v].second;
    tree.cut(v, neighbor1);
    tree.cut(v, neighbor2);
    chain_map.erase(v);
    if(chain_map[neighbor1].first == v) {chain_map[neighbor1].first = neighbor2;} else {chain_map[neighbor1].second = neighbor2;}
    if(chain_map[neighbor2].first == v){chain_map[neighbor2].first = neighbor1;} else {chain_map[neighbor2].second = neighbor1;}
    tree.link(neighbor1, neighbor2, id);
  }
  head_map[v - n] = MAX_VERTEX_T;
}

/**********************/

template<typename DynamicTree, typename aug_t>
void TernarizedTree<DynamicTree, aug_t>::link(vertex_t u, vertex_t v, aug_t weight) {
  // If vertex degree <= 3, proceed with normal link, otherwise ternarize vertex
  // and add new vertex to head.

  //TODO: Change get_degree method of RCTree to have default value for round=0
  //TODO: Add a get_neighbors method to both Topology and RC trees.
  
  assert(!connected(u,v));
  if(u > v) std::swap(u,v);
  if(tree.get_degree(u) < 3 && tree.get_degree(v) < 3) {
    tree.link(u,v, weight); 
    edge_map[std::pair<vertex_t, vertex_t>(u,v)] = std::pair<vertex_t, vertex_t>(u,v);
    return;
  }
  vertex_t to_link_u = u;
  vertex_t to_link_v = v;
  if(tree.get_degree(u) == 3) to_link_u = determine_link_v(u); 
  if(tree.get_degree(v) == 3) to_link_v = determine_link_v(v); 
  tree.link(to_link_u, to_link_v, weight);
  edge_map[std::pair<vertex_t, vertex_t>(u,v)] = std::pair<vertex_t, vertex_t>(to_link_u, to_link_v);
}

template<typename DynamicTree, typename aug_t>
void TernarizedTree<DynamicTree, aug_t>::cut(vertex_t u, vertex_t v) {
  // Find vertex on ternarized path that points to v and u
  // Remove this from the ternarized path.
  if(u > v){std::swap(u,v);}

  // Test of hash function to make sure we are able to find pairs in the map.
  if(edge_map.find(std::pair(u,v)) == edge_map.end()){
    std::cerr << "(u,v) not in map\n";
    pair_hash p;
    std::cerr << p(std::pair(u,v)) << "\n";
    for(auto pair : edge_map){
      std::cerr << pair.first.first << "," << pair.first.second << " : " << pair.second.first << "," << pair.second.second << " ";
      std::cerr << p(pair.first) << " " << (p(pair.first) == p(std::pair(u,v))) << "\n";
    }
    assert(false);
  }
  auto edge_pair = edge_map[std::pair(u,v)];
  vertex_t v1 = edge_pair.first, v2 = edge_pair.second;
  tree.cut(v1, v2);
  if(v1 >= n){
    delete_ternarized_vertex(v1);
    free_ids.push(v1);
  }
  if(v2 >= n){
    delete_ternarized_vertex(v2);
    free_ids.push(v2);
  }
  edge_map.erase(std::pair(u,v));
}

template<typename DynamicTree, typename aug_t>
bool TernarizedTree<DynamicTree, aug_t>::connected(vertex_t u, vertex_t v) {
  return tree.connected(u,v);
}

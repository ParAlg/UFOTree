#include <unordered_map>
#include <queue>
#include "types.h"
#include <cassert>
#include "util.h"


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
private:
  vertex_t determine_link_v(vertex_t u);
  void ternarize_vertex(vertex_t u);
  void delete_ternarized_vertex(vertex_t v);
  vertex_t get_id(vertex_t u);
  vertex_t get_edge_val(vertex_t u);
  // Underlying dynamic tree data structure
  DynamicTree tree;
  // Ternarization book-keeping
  vertex_t n;
  aug_t id;
  std::unordered_map<std::pair<vertex_t, vertex_t>, std::pair<vertex_t, vertex_t>> edge_map;
  std::queue<vertex_t> free_ids;
};

template<typename DynamicTree, typename aug_t>
TernarizedTree<DynamicTree, aug_t>::TernarizedTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d) : 
  tree(2*n, q, f, id, d), n(n) {}


/*** HELPER METHODS ***/

template<typename DynamicTree, typename aug_t>
vertex_t TernarizedTree<DynamicTree, aug_t>::determine_link_v(vertex_t u){
  auto u_neighbors = tree.get_neighbors(u);
  // Find ternarized neighbor if any exists.
  auto ternary_neighbor_u = nullptr; 
  for(auto neighbor_c : u_neighbors){if(get_id(neighbor_c) >= n) ternary_neighbor_u = neighbor_c; break;}
  if(ternary_neighbor_u == nullptr){
    ternarize_vertex(u); 
    return u;
  } else {
    tree.cut(u, get_id(ternary_neighbor_u));
    vertex_t new_ternarized_n = free_ids.front(); free_ids.pop();
    tree.link(u, new_ternarized_n, id);
    return new_ternarized_n;
  }
}

template<typename DynamicTree, typename aug_t>
void TernarizedTree<DynamicTree, aug_t>::ternarize_vertex(vertex_t u){
  auto u_neighbors = tree.get_neighbors(u);
  auto to_cut = get_id(u_neighbors[0]);
  auto edge_val = get_edge_val(u);
  tree.cut(u, to_cut);
  vertex_t new_ternarized_n = free_ids.front(); free_ids.pop();
  tree.link(u, new_ternarized_n, id);
  tree.link(new_ternarized_n, to_cut, edge_val);
}

template<typename DynamicTree, typename aug_t>
void TernarizedTree<DynamicTree, aug_t>::delete_ternarized_vertex(vertex_t v){
  if(tree.get_degree(v) == 1){
    vertex_t neighbor = MAX_VERTEX_T;
    for(int i = 0; i < 3; i++){if(tree.get_neighbors(v)[i] != nullptr) get_id(tree.get_neighbors(v)[i]);}
    tree.cut(v, neighbor);
  } else {
    vertex_t n1 = MAX_VERTEX_T, n2 = MAX_VERTEX_T;
    for(int i = 0; i < 3; i++){
      if(tree.get_neighbors[v][i] != nullptr){
        // 1 liner or 5 liner?
        if(n1 == MAX_VERTEX_T) {n1 = get_id(tree.get_neighbors[v][i]);} else {n2 = get_id(tree.get_neighbors[v][i]);}
      }
    }
    tree.cut(v, n1);
    tree.cut(v, n2);
    tree.link(n1, n2, id);
  }
}

template<typename DynamicTree, typename aug_t>
vertex_t TernarizedTree<DynamicTree, aug_t>::get_id(vertex_t u){}

template<typename DynamicTree, typename aug_t>
vertex_t TernarizedTree<DynamicTree, aug_t>::get_edge_val(vertex_t u){}
/**********************/
template<typename DynamicTree, typename aug_t>
void TernarizedTree<DynamicTree, aug_t>::link(vertex_t u, vertex_t v, aug_t weight) {
  // If vertex degree <= 3, proceed with normal link, otherwise ternarize vertex
  // and add new vertex to head.

  //TODO: Change get_degree method of RCTree to have default value for round=0
  //TODO: Add a get_neighbors method to both Topology and RC trees.

  if(tree.get_degree(u) < 3 && tree.get_degree(v) < 3) 
    tree.link(u,v, weight); return;

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
  assert(((void)"(u,v) not in map", edge_map.find(std::pair(u,v)) != edge_map.end()));
  auto edge_pair = edge_map.find(std::pair(u,v))->second;
  vertex_t v1 = edge_pair.first, v2 = edge_pair.second; 
  tree.cut(v1, v2);

  if(get_id(v1) >= n){
    delete_ternarized_vertex(v1);
    free_ids.push(v1);
  }

  if(get_id(v2) >= n){
    delete_ternarized_vertex(v2);
    free_ids.push(v2);
  }
}

template<typename DynamicTree, typename aug_t>
bool TernarizedTree<DynamicTree, aug_t>::connected(vertex_t u, vertex_t v) {
  return tree.connected(u,v);
}

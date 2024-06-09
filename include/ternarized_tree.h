#include <vector>
#include "types.h"
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
    // Underlying dynamic tree data structure
    DynamicTree tree;
    // Ternarization book-keeping
    vertex_t n;
    std::vector<vertex_t> free_vertex_ids;
    // Store the head of the ternarized vertices linked list for each ternarized vertex.
    std::vector<vertex_t> vertex_head_map;  
};

template<typename DynamicTree, typename aug_t>
TernarizedTree<DynamicTree, aug_t>::TernarizedTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d) : 
    tree(2*n, q, f, id, d), n(n) {}

template<typename DynamicTree, typename aug_t>
void TernarizedTree<DynamicTree, aug_t>::link(vertex_t u, vertex_t v, aug_t weight) {
  // If vertex degree <= 3, proceed with normal link, otherwise ternarize vertex
  // and add new vertex to head.
  
  //TODO: Change get_degree method of RCTree to have default value for round=0
  //TODO: Add a get_neighbors method to both Topology and RC trees.

  if(tree.get_degree(u) < 3 && tree.get_degree(v) < 3) 
    tree.link(u,v, weight); return;

  if(tree.get_degree(u) == 3){
    auto u_neighbors = tree.get_neighbors(u);
    auto ternary_neighbor = nullptr;
    for(auto neighbor_c : u_neighbors){
      auto neighbor = neighbor_c->v;
      if(neighbor> n) ternary_neighbor = neighbor_c; break;
    }
    
  }
}

template<typename DynamicTree, typename aug_t>
void TernarizedTree<DynamicTree, aug_t>::cut(vertex_t u, vertex_t v) {
  // Find vertex on ternarized path that points to v and u
  // Remove this from the ternarized path. 
  
    tree.cut(u, v);
}

template<typename DynamicTree, typename aug_t>
bool TernarizedTree<DynamicTree, aug_t>::connected(vertex_t u, vertex_t v) {
    return tree.connected(u,v);
}

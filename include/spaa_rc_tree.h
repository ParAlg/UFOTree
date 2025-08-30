#pragma once
#include "../spaa_rc_tree/RCtrees/RC.h"
#include "../spaa_rc_tree/RCtrees/RCdynamic.h"
#include "../spaa_rc_tree/RCtrees/ternarizer.h"
#include "types.h"
#include "../spaa_rc_tree/RCtrees/RC_test.h"
#include "../spaa_rc_tree/RCtrees/path_query.h"
#include <limits>
#include <tuple>

using namespace dgbs;


size_t static_space_used = 0;

template<typename aug_t>
class ParallelRCTree{
  // Summarize how this code is organized?
  // RC tree code in RC.h
  // Insertion code in RCdynamic.h
  // Query code in path_query and subtree_query
  //
  // Ternarizer stores a simpler smaller version of the original tree to determine what nodes to ternarize
  // Insertion edges are first passed into the ternarizer which mends then outputs a list of ternarized insertion edges
  // these are then inserted into the actual RC tree
  //
  // Space - Things in RC, RCDynamic, and ternarizer file.
 
public:
  parlay::sequence<cluster<int,aug_t> > clusters;
  ternarizer<int, int> tr = ternarizer(200, std::numeric_limits<int>::min()); // NOTE: tern_node uses -1 as default edge array index - need to fix that before converting this to int
  int n;
  int k; 
  void batch_link(parlay::sequence<std::tuple<int, int, aug_t>>& links);
  void batch_cut(parlay::sequence<std::pair<int,int>> & cuts);
  void link(int u, int v, aug_t w=std::numeric_limits<aug_t>::min());
  void cut(int u, int v);
  size_t space();
  aug_t identity;
  aug_t path_query(int u, int v);
  std::function<aug_t(aug_t,aug_t)> func;
  void verify_tree_correctness(); 

  ParallelRCTree(int _n, int _k = 1, aug_t id = std::numeric_limits<int>::min(), std::function<aug_t(aug_t,aug_t)> _func = [] (int A, int B) {return std::max(A,B);}){
    n = _n;
    k = _k;
    identity = id;
    parlay::sequence<std::tuple<int,int, int> > initial_edges;
    create_base_clusters(clusters, initial_edges, static_cast<int>(3), n);
    create_RC_tree(clusters, n, 0, _func); 
    func = _func;
  }
  ~ParallelRCTree(){
    deleteRCtree(clusters);
    static_space_used = 0;
  }
};

template<typename aug_t>
void ParallelRCTree<aug_t>::link(int u, int v, aug_t w){
  w = identity;
  parlay::sequence<std::tuple<int,int,aug_t> > insert_edges; 
  insert_edges.push_back({u,v,w});
  parlay::sequence<std::pair<int,int>> delete_edges;
  batchInsertEdge(delete_edges, insert_edges, clusters, 0, func); 
}

template<typename aug_t>
void ParallelRCTree<aug_t>::cut(int u, int v){
  parlay::sequence<std::tuple<int,int,aug_t> > insert_edges; 
  parlay::sequence<std::pair<int,int>> delete_edges; delete_edges.push_back({u,v});
  batchInsertEdge(delete_edges, insert_edges, clusters, 0, func); 
}
template<typename aug_t>
void ParallelRCTree<aug_t>::batch_link(parlay::sequence<std::tuple<int,int, aug_t>>& links){ 
  parlay::sequence<std::pair<int,int>> delete_edges;
  batchInsertEdge(delete_edges, links, clusters, 0, func); 
}

template<typename aug_t>
void ParallelRCTree<aug_t>::batch_cut(parlay::sequence<std::pair<int,int>>& cuts){
  parlay::sequence<std::tuple<int, int, aug_t>> add_edges;
  batchInsertEdge(cuts, add_edges, clusters, 0, func); 
}

template<typename aug_t>
aug_t ParallelRCTree<aug_t>::path_query(int u, int v){
  return pathQuery<int,aug_t, std::function<aug_t(aug_t,aug_t)> >(u,v,clusters,identity, func) ;
}

template<typename aug_t>
void ParallelRCTree<aug_t>::verify_tree_correctness(){
  check_children_values(clusters);
  check_parents_children(clusters);
}

template<typename aug_t>
size_t ParallelRCTree<aug_t>::space(){
  auto ans = parlay::type_allocator<cluster<int,aug_t>>::num_used_bytes() + parlay::type_allocator<node<int,aug_t>>::num_used_bytes() + static_space_used;
  return ans; 
}
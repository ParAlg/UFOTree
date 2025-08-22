#pragma once
#include "../spaa_rc_tree/RCtrees/RC.h"
#include "../spaa_rc_tree/RCtrees/RCdynamic.h"
#include "../spaa_rc_tree/RCtrees/ternarizer.h"
#include "types.h"
#include "../spaa_rc_tree/RCtrees/RC_test.h"
#include <limits>
#include <tuple>

using namespace dgbs;

template<typename aug_t>
class ParallelRCTree{
  // Summarize how this code is organized?
  // RC tree code in RC.h
  // Insertion code in RCdynamic.h
  // Query code in path_query and subtree_query
  //
  // Ternarizer stores a simpler smaller version of the original tree to determine what nodes to ternarize
  // Insertion/Deletion edges are first passed into the ternarizer which mends then outputs a list of ternarized and deletion edges (i.e. 
  // what edges need to be inserted/deleted in the underlying ternarized tree)
  // these are then inserted into the actual RC tree via the batchInsertEdge method 
  //
  // Space - Things in RC, RCDynamic, and ternarizer file.
 
public:
  parlay::sequence<cluster<int,aug_t> > clusters;
  //ternarizer<int, int> tr = ternarizer(200, std::numeric_limits<int>::min()); // NOTE: tern_node uses -1 as default edge array index - need to fix that before converting this to int
  int n;
  int k; 
  void batch_link(parlay::sequence<std::tuple<int, int, aug_t>>& links);
  void batch_cut(parlay::sequence<std::pair<int,int>> & cuts);
  aug_t path_query(int u, int v);
  std::function<aug_t(aug_t,aug_t)> func;
  void verify_tree_correctness();
  ParallelRCTree(int _n, int _k = 1, std::function<aug_t(aug_t,aug_t)> _func = [] (int A, int B) {return std::max(A,B);}){
    n = _n;
    k = _k;
    //tr = ternarizer(_n * extra_tern_node_factor, std::numeric_limits<int>::min());
    parlay::sequence<std::tuple<int,int, int>> initial_edges;
    create_base_clusters(clusters, initial_edges, static_cast<int>(3), n);
    create_RC_tree(clusters, n, 0, _func); 
    func = _func;
  }
};

template<typename aug_t>
void ParallelRCTree<aug_t>::batch_link(parlay::sequence<std::tuple<int,int, aug_t>>& links){ 
  parlay::sequence<std::pair<int,int>> delete_edges; // Created for a valid call to batchInsertEdge, does not modify the tree at all.
  //auto ternarizedEdges = tr.add_edges(links);
  batchInsertEdge(delete_edges, links, clusters, 0, func); 
}

template<typename aug_t>
void ParallelRCTree<aug_t>::batch_cut(parlay::sequence<std::pair<int,int>>& cuts){
  parlay::sequence<std::tuple<int, int, aug_t>> add_edges;
  batchInsertEdge(cuts, add_edges, clusters, 0, func); 
}

template<typename aug_t>
aug_t ParallelRCTree<aug_t>::path_query(int u, int v){
  
}

template<typename aug_t>
void ParallelRCTree<aug_t>::verify_tree_correctness(){
  check_children_values(clusters);
  check_parents_children(clusters);
}


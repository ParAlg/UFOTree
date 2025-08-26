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
class ParallelRCTreeTernarized{
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
  ternarizer<int,int>* tr;// NOTE: tern_node uses -1 as default edge array index - need to fix that before converting this to int
  int n;
  int k; 
  aug_t identity;
  void batch_link(parlay::sequence<std::tuple<int, int, aug_t>>& links);
  void batch_cut(parlay::sequence<std::pair<int,int>> & cuts);
  void link(int u, int v, aug_t w=std::numeric_limits<aug_t>::min());
  void cut(int u, int v);
  aug_t path_query(int u, int v);
  std::function<aug_t(aug_t,aug_t)> func;
  void verify_tree_correctness();
  ParallelRCTreeTernarized(int _n, int _k = 1, aug_t id = std::numeric_limits<int>::min(), std::function<aug_t(aug_t,aug_t)> _func = [] (int A, int B) {return std::max(A,B);}){
    n = _n * extra_tern_node_factor; 
    //cout << "N: " << n << "\n";
    k = _k;
    identity = id;
    tr = new ternarizer<int,int>{_n, id};
    parlay::sequence<std::tuple<int, int, int> > initial_edges;
    create_base_clusters(clusters, initial_edges, static_cast<int>(3), n);
    create_RC_tree(clusters, n, 0, _func); 
    func = _func;
  }
  ~ParallelRCTreeTernarized(){
    free(tr);
    deleteRCtree(clusters);
  }
};

template<typename aug_t>
void ParallelRCTreeTernarized<aug_t>::link(int u, int v, aug_t w){
  w = identity;
  parlay::sequence<std::tuple<int,int,aug_t> > insert_edges; 
  insert_edges.push_back({u,v,w});
  batch_link(insert_edges);
}

template<typename aug_t>
void ParallelRCTreeTernarized<aug_t>::cut(int u, int v){
  parlay::sequence<std::pair<int,int>> cuts; 
  cuts.push_back({u,v});
  batch_cut(cuts);
}

template<typename aug_t>
void ParallelRCTreeTernarized<aug_t>::batch_link(parlay::sequence<std::tuple<int,int, aug_t>>& links){ 
  auto ternarizedEdges = tr->add_edges(links);
  /*std::cout << "Delete Edges :" << "\n";
  for(auto e : ternarizedEdges.first){
    std::cout << e.first << " " << e.second << "\n";
  }
  std::cout << "Insert Edges :" << "\n";*/
  /*for(auto &e : ternarizedEdges.first){
    if(std::get<0>(e) >= n){
      std::get<0>(e) -= n/2;
    }
    if(std::get<1>(e) >= n){
      std::get<1>(e) -= n/2;
    }
    //std::cout << std::get<0>(e) << " " << std::get<1>(e) << "\n"; 
  }

  for(auto &e : ternarizedEdges.second){
    if(std::get<0>(e) >= n){
      std::get<0>(e) -= n/2;
    }
    if(std::get<1>(e) >= n){
      std::get<1>(e) -= n/2;
    }
    //std::cout << std::get<0>(e) << " " << std::get<1>(e) << "\n"; 
  }
  /*for(auto e : ternarizedEdges.second){
    std::cout << std::get<0>(e) << " " << std::get<1>(e) << "\n"; 
  }*/
  //std::cout << "\n";
  //parlay::sequence<std::tuple<int,int,int> > s; s.push_back(ternarizedEdges.second[0]);*/
  batchInsertEdge(ternarizedEdges.first, ternarizedEdges.second, clusters, 0, func); 
}

template<typename aug_t>
void ParallelRCTreeTernarized<aug_t>::batch_cut(parlay::sequence<std::pair<int,int>>& cuts){
  auto ternarizedEdges = tr->delete_edges(cuts);
  /*for(auto &e : ternarizedEdges.first){
    if(std::get<0>(e) >= n){
      std::get<0>(e) -= n/2;
    }
    if(std::get<1>(e) >= n){
      std::get<1>(e) -= n/2;
    }
    //std::cout << std::get<0>(e) << " " << std::get<1>(e) << "\n"; 
  }
  for(auto &e : ternarizedEdges.second){
    if(std::get<0>(e) >= n){
      std::get<0>(e) -= n/2;
    }
    if(std::get<1>(e) >= n){
      std::get<1>(e) -= n/2;
    }
    //std::cout << std::get<0>(e) << " " << std::get<1>(e) << "\n"; 
  }*/
  batchInsertEdge(ternarizedEdges.first, ternarizedEdges.second, clusters, 0, func); 
}

template<typename aug_t>
aug_t ParallelRCTreeTernarized<aug_t>::path_query(int u, int v){
  return pathQuery<int,aug_t, std::function<aug_t(aug_t,aug_t)> >(u,v,clusters,identity, func);
}

template<typename aug_t>
void ParallelRCTreeTernarized<aug_t>::verify_tree_correctness(){
  check_children_values(clusters);
  check_parents_children(clusters);
}


#pragma once
#include "types.h"
#include "util.h"
#include "ufo_cluster.h"
#include <absl/container/flat_hash_set.h>
#include <absl/container/flat_hash_map.h>
/*
 * */

template<typename v_t, typename e_t>
class IUFOCluster;

template <typename v_t, typename e_t>
class IUFOCluster{
public:
  // Data
  v_t id;  
  int level;
  v_t parent;
  IUFOCluster<v_t,e_t>* next, *prev;  // Fields for the linked list in the child data structure.
  v_t tail; 
  std::tuple<v_t,v_t,e_t> contracting_edge, renamed_edge;
  // Methods 

  // Disconnect all children above level l
  void add_child(IUFOCluster* child);
  IUFOCluster(int _id=-1, int _level=0, v_t _parent=-1, 
               std::tuple<v_t,v_t,e_t> _contracting_edge= std::make_tuple(v_t(), v_t(),e_t()), 
               std::tuple<v_t,v_t,e_t> _renamed_edge= std::make_tuple(v_t(), v_t(), e_t()),
               IUFOCluster<v_t, e_t> *_next = nullptr, IUFOCluster<v_t, e_t> *_prev = nullptr,
               v_t _tail = -1){
    id = _id;
    level = _level;
    parent = _parent;
    contracting_edge = _contracting_edge;
    renamed_edge = _renamed_edge;
    next = _next;
    prev = _prev;
    tail = _tail;
  }
};

template<typename v_t, typename e_t>
class IUFOTree{
// Made Public for testing purposes
public:
  // Data
  int n;
  std::vector<IUFOCluster<v_t, e_t>*> clusters;
  absl::flat_hash_map<std::tuple<v_t,v_t,e_t>, std::vector<std::tuple<v_t,v_t,e_t> > > edge_history;
  std::vector<int> high_deg_level;
  std::vector<std::vector<v_t> > root_clusters;

  // Methods  

  void link(v_t u, v_t v, e_t w);
  void cut(v_t u, v_t v);
  void recluster_tree();
  e_t path_query(v_t u, v_t v);
  bool is_valid();
  void print_tree();
  void remove_ancestors(v_t v, v_t u);   
  void disconnect_children(IUFOCluster<v_t, e_t> c);
  IUFOTree (int _n = 0){
    n = _n;
    high_deg_level.resize(n);
    clusters.resize(n);
    root_clusters.resize(100);
  }
};

template<typename v_t, typename e_t>
void IUFOTree<v_t, e_t>::disconnect_children(IUFOCluster<v_t, e_t> c){
  if (c.tail != -1){
    auto curr = clusters[c.tail];
    while(curr != nullptr && curr->level >= c.level){
      curr->parent = -1;
      root_clusters[curr->level].push_back(curr->id);
      curr->next = nullptr;
      auto prev = curr->prev;
      curr->prev = nullptr;
      curr = prev;
    }
    if(curr != nullptr){
      c.tail = curr->id;
    }else{
      c.tail = -1;
    }
  }
}

template<typename v_t, typename e_t>
void IUFOTree<v_t, e_t>::remove_ancestors(v_t u, v_t v){ 
  auto cluster = *clusters[u];
  int prev_level = 0;
  while(cluster.parent != -1){
    auto tmp = cluster.level;
    cluster.level = std::max(prev_level, high_deg_level[cluster.id]); //NOTE: SHOULD WE DEREFERENCE HERE OR STORE EVERYTHING IN MEMORY 
    prev_level = tmp;
    disconnect_children(cluster); 
    auto prev_parent = *clusters[cluster.parent];
    cluster.parent = -1;
    cluster = prev_parent;
  }

  cluster = *clusters[v];
  prev_level = 0;
  while(cluster.parent != -1){ 
    auto tmp = cluster.level;
    cluster.level = std::max(prev_level, high_deg_level[cluster.id]); //NOTE: SHOULD WE DEREFERENCE HERE OR STORE EVERYTHING IN MEMORY 
    prev_level = tmp;
    disconnect_children(cluster);
    auto prev_parent = *clusters[cluster.parent];
    cluster.parent = -1;
    cluster = prev_parent;
  }
  //For the final cluster on the ancestor path of v
  cluster.level = std::max(prev_level, high_deg_level[cluster.id]); 
  disconnect_children(cluster);
}

template<typename v_t, typename e_t>
void IUFOTree<v_t, e_t>::recluster_tree(){
  
}

template<typename v_t, typename e_t>
void IUFOTree<v_t, e_t>::link(v_t u, v_t v, e_t w){
  remove_ancestors(u);   
  remove_ancestors(v);
  recluster_tree();
}

template<typename v_t, typename e_t>
void IUFOTree<v_t, e_t>::cut(v_t u, v_t v){

}

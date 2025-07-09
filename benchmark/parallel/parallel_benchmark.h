#pragma once
#include <vector>
#include <iostream>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include <parlay/internal/get_time.h>
#include "types.h"
#include <spaa_rc_tree.h>

namespace parallel_dynamic_tree_benchmark {

// Returns the time in seconds to perform all of the updates
template <typename DynamicTree>
double get_update_speed(vertex_t n, vertex_t k, std::vector<std::vector<Update>> update_sequences) {
  parlay::internal::timer my_timer("");
  DynamicTree t(0,0);
  ParallelRCTree<int> type_finder(0,0);
  if(typeid(t).name() == typeid(type_finder).name()){  
    for(auto updates : update_sequences){  
      DynamicTree tree(n,k);
      parlay::sequence<std::pair<vertex_t, vertex_t> > deletion_batch;
      parlay::sequence<std::tuple<vertex_t, vertex_t, int> > insertion_batch;
      UpdateType batch_type = updates[0].type; 
      my_timer.start();
      for (auto update : updates) {
        if (update.type != batch_type) {
          if(batch_type==INSERT) {
            tree.batch_link(insertion_batch);
            insertion_batch.clear();
          } else {
            tree.batch_cut(deletion_batch);
            deletion_batch.clear();
          }
          batch_type = updates[0].type;
        }
        if(batch_type == INSERT){
          insertion_batch.push_back(std::make_tuple(update.edge.src, update.edge.dst, 0)); 
        } else{
          deletion_batch.push_back(std::make_pair(update.edge.src, update.edge.dst));
        }
        if(insertion_batch.size() == k || deletion_batch.size() == k) {
          if(batch_type==INSERT) {
            tree.batch_link(insertion_batch);
            insertion_batch.clear();
          } else {
            tree.batch_cut(deletion_batch);
            deletion_batch.clear();
          }
        }
      }
      if(batch_type==INSERT) {
            tree.batch_link(insertion_batch);
            insertion_batch.clear();
          } else {
            tree.batch_cut(deletion_batch);
            deletion_batch.clear();
          }
      my_timer.stop();
    }
  } else{
    for (auto updates : update_sequences) {
      DynamicTree tree(n, k);
      parlay::sequence<Edge> batch;
      UpdateType batch_type = updates[0].type;
      my_timer.start();
      for (auto update : updates) {
        if (update.type != batch_type) {
          (batch_type==INSERT) ? tree.batch_link(batch) : tree.batch_cut(batch);
          batch.clear();
          batch_type = updates[0].type;
        }
        batch.push_back(update.edge);
        if (batch.size() == k) {
          (batch_type==INSERT) ? tree.batch_link(batch) : tree.batch_cut(batch);
          batch.clear();
        }
      }
      (batch_type==INSERT) ? tree.batch_link(batch) : tree.batch_cut(batch);
      my_timer.stop();
    }
  }
  return my_timer.total_time()/update_sequences.size();
}

}

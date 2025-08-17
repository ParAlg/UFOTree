#pragma once
#include "parlay/parallel.h"
#include "types.h"
#include "util.h"
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <parlay/internal/get_time.h>
#include <random>
#include <utility>


/*Benchmarking:
- Batches: 1, 10, 100, 1000, 10^6, 10^7
- Measure across different cores
- Randomized insertions, deletions - situations where we can’t really take advantage of cache coherency
*/

using namespace dgbs;
uint32_t batch_size = 1000;
//uint32_t single_size = 100;
uint32_t granularity = 100;

// Try out random insertions and deletions vs batched insertions/deletions

template <typename par_cluster>
static void benchmark_single_insert_delete_randomized(){
  parlay::internal::timer my_timer("");
  par_cluster cluster(granularity);
  parlay::sequence<std::pair<par_cluster*, bool> > updates;
  for(int i = 0; i < batch_size; i++){
    par_cluster new_cluster;
    cluster.insert_neighbor(&new_cluster);
    updates.push_back(std::make_pair(&new_cluster, 1));
    par_cluster insert_cluster;
    updates.push_back(std::make_pair(&insert_cluster, 0));
  }
  std::shuffle(updates.begin(), updates.end(), std::mt19937());
  
  my_timer.start();
  parlay::parallel_for(0, 2*batch_size, [&] (int i){
    if(updates[i].second){
      cluster.delete_neighbor(updates[i].first);
    } else{
      cluster.insert_neighbor(updates[i].first);
    }
  }, granularity);
  my_timer.stop();
  std::cout << "Total Time for " << (2 * batch_size) << " Inserts and Deletes: " << my_timer.total_time() << "\n";
}

template <typename par_cluster>
static void benchmark_batch_updates(){
  parlay::internal::timer my_timer("");
  parlay::sequence<par_cluster*> clusters;
  
  for(int i = 1; i <= batch_size; i++){
    par_cluster* new_cluster = new par_cluster;
    clusters.push_back(new_cluster);
  }
  par_cluster cluster(granularity);
  my_timer.start();
  cluster.insert_neighbors(clusters);
  my_timer.stop();
  auto insertion_time = my_timer.total_time(); 
  my_timer.reset();
  
  my_timer.start();
  cluster.delete_neighbors(clusters);
  my_timer.stop();
  auto deletion_time = my_timer.total_time();

  std::cout << "Batch Insertion Time: " << insertion_time << "\n";
  std::cout << "Batch Deletion Time: " << deletion_time << "\n";
}

template<typename par_cluster>
static void benchmark_real_trees(int n, std::function<std::vector<Update>(vertex_t, long)> update_fn,
                                  bool serial=true, bool batch=true){
  
  auto update_sequences = update_fn(n, rand());
  std::vector<par_cluster> clusters(n);
  std::vector<Update> insertions, deletions;
  int delete_idx = -1;
  for(int i = 0; i < update_sequences.size(); i++){
    if(update_sequences[i].type == INSERT){
      insertions.push_back(update_sequences[i]);
    } else{
      deletions.push_back(update_sequences[i]);
    }
  }
  
  parlay::internal::timer my_timer;
  if(serial){
    // Serial Insertions
    my_timer.start();
    parlay::parallel_for(0, insertions.size(), [&] (int i){
      auto edge = insertions[i].edge;
      clusters[edge.src].insert_neighbor(&clusters[edge.dst]);
    }, granularity);
    my_timer.stop();
    auto insertion_time = my_timer.total_time();
    my_timer.reset();
    my_timer.start();
    parlay::parallel_for(0, deletions.size(), [&] (int i ){
      auto edge = deletions[i].edge;
      clusters[edge.src].delete_neighbor(&clusters[edge.dst]);
    }, granularity);
    my_timer.stop();  
    auto deletion_time = my_timer.total_time();
    std::cout << "Serial Insertion Time:" << insertion_time << "\n";
    std::cout << "Serial Deletion Time: " << deletion_time << "\n";
    my_timer.reset();
  }
  if(batch){
    my_timer.start();
    parlay::sort_inplace(insertions, [&](Update a, Update b){
      return a.edge.src < b.edge.src;
    });

    auto starts = parlay::delayed_tabulate(insertions.size(), [&] (int i) {
        if (i == 0 || insertions[i].edge.src != insertions[i-1].edge.src) return true;
        return false;
    });
    auto offsets = parlay::pack_index(starts);
    auto final_updates = parlay::tabulate(offsets.size(), [&] (size_t i) {
        size_t start = offsets[i];
        size_t end = i == offsets.size()-1 ? insertions.size() : offsets[i+1];
        return std::make_pair(insertions[offsets[i]].edge.src, parlay::make_slice(insertions.begin() + start, insertions.begin() + end));
    });

    parlay::parallel_for(0, final_updates.size(), [&] (size_t i) {
            auto& [cluster, edges] = final_updates[i];
            auto neighbors = parlay::map(edges, [&] (auto x) { return &clusters[x.edge.dst]; });
            clusters[i].insert_neighbors(neighbors);
    });
    my_timer.stop();
    std::cout << "Parallel Insertion: " << my_timer.total_time() << "\n";
    my_timer.reset();

    my_timer.start();
    parlay::sort_inplace(deletions, [&](Update a, Update b)
                         { return a.edge.src < b.edge.src; });

    auto new_starts = parlay::delayed_tabulate(deletions.size(), [&](int i){
        if (i == 0 || deletions[i].edge.src != deletions[i-1].edge.src) return true;
        return false; });
    offsets = parlay::pack_index(new_starts);
    final_updates = parlay::tabulate(offsets.size(), [&](size_t i){
        size_t start = offsets[i];
        size_t end = i == offsets.size()-1 ? deletions.size() : offsets[i+1];
        return std::make_pair(deletions[offsets[i]].edge.src, parlay::make_slice(deletions.begin() + start, deletions.begin() + end)); });

    parlay::parallel_for(0, final_updates.size(), [&](size_t i){
            auto& [cluster, edges] = final_updates[i];
            auto neighbors = parlay::map(edges, [&] (auto x) { return &clusters[x.edge.dst]; });
            clusters[i].delete_neighbors(neighbors); });
    my_timer.stop();
    std::cout << "Parallel Deletion: " << my_timer.total_time() << "\n";
  }
}
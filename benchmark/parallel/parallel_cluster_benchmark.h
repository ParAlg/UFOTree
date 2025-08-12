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

// This is probably useless:

/*
template <typename par_cluster>
static void benchmark_single_insert_delete_continuous(){
  parlay::internal::timer my_timer("");
  parlay::sequence<par_cluster*> insertions;
  insertions.resize(100);
  par_cluster cluster(granularity);
  my_timer.start();
  parlay::parallel_for(0, single_size, [&](int i){
    cluster.insert_neighbor(insertions[i]);
  });
  my_timer.stop();
  auto insertion_time = my_timer.total_time();
  my_timer.reset();

  std::vector<int> deletions;
  for(int i = 0; i < single_size; i++){
    deletions.push_back(i);
  }
  std::shuffle(deletions.begin(), deletions.end(), std::mt19937());
  my_timer.start();
  parlay::parallel_for(0, single_size, [&](int i){
    cluster.delete_neighbor(insertions[i]);
  });
  my_timer.stop();
  auto deletion_time = my_timer.total_time();

  std::cout << "Insertion Time: " << insertion_time << "\n" << "Deletion Time: " << deletion_time;
}
*/


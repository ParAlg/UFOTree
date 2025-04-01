#pragma once
#include <vector>
#include <iostream>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include <parlay/internal/get_time.h>
#include "types.h"
#include "graph_utils.h"
#include "util.h"


namespace graph_benchmark {

std::vector<Update> stream_file_BFS_benchmark(parlay::sequence<parlay::sequence<vertex_t>>& G, long seed = -1) {
    if (seed == -1) seed = time(NULL);
    srand(seed);
    auto E = graph_utils::BFS_forest(G, seed);
    auto insertions = parlay::map(parlay::random_shuffle(E, parlay::random(rand())), [&](auto e) -> Update {return {INSERT, e};});
    auto deletions = parlay::map(parlay::random_shuffle(E, parlay::random(rand())), [&](auto e) -> Update {return {DELETE, e};});
    return parlay::append(insertions, deletions).to_vector();
}

std::vector<Update> stream_file_RIS_benchmark(parlay::sequence<parlay::sequence<vertex_t>>& G, long seed = -1) {
    if (seed == -1) seed = time(NULL);
    srand(seed);
    auto E = graph_utils::RIS_forest(G, seed);
    auto insertions = parlay::map(parlay::random_shuffle(E, parlay::random(rand())), [&](auto e) -> Update {return {INSERT, e};});
    auto deletions = parlay::map(parlay::random_shuffle(E, parlay::random(rand())), [&](auto e) -> Update {return {DELETE, e};});
    return parlay::append(insertions, deletions).to_vector();
}

template<typename DynamicTree>
double incremental_MSF_benchmark(parlay::sequence<std::pair<int, Edge>>& edges, int n) {
  Edge garbage; garbage.src = MAX_VERTEX_T; garbage.dst = MAX_VERTEX_T;
  std::pair<int, Edge> id(std::numeric_limits<int>::min(), garbage);
  DynamicTree tree(n, PATH, [] (std::pair<int, Edge> a, std::pair<int, Edge> b){return a.first > b.first ? a : b;}, id, id);
  parlay::internal::timer my_timer("");
  
  my_timer.start();
  for (auto edge : edges) {
    vertex_t u = edge.second.src;
    vertex_t v = edge.second.dst;
    if (!tree.connected(u, v)) {
      tree.link(u, v, edge);
      continue;
    }
    std::pair<int, Edge> max_weight_edge = tree.path_query(u, v);
    if (max_weight_edge.first > edge.first) {
      tree.cut(max_weight_edge.second.src, max_weight_edge.second.dst);
      tree.link(u, v, edge);
    }
  }

  my_timer.stop();
  return my_timer.total_time();
}

double incremental_conn_benchmark(parlay::sequence<parlay::sequence<vertex_t>>& G, long seed = -1){
  return 0.0;
}
};

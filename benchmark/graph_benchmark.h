#pragma once
#include <vector>
#include <iostream>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include <parlay/internal/get_time.h>
#include "types.h"
#include "graph_utils.h"


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

/*
template<typename DynamicTree>
double incremental_MSF_benchmark(std::vector<std::pair<int, Edge>>& edges, int n){
  DynamicTree tree(n, PATH, [] (std::pair<int, Edge> a, std::pair<int, Edge> b){return a.first > b.first ? a : b;},  std::numeric_limits<int>::max(), 0);

  for(edge_pair : (*edges)){
    Edge e = edge_pair.first;
    int weight = edge_pair.second;
    
    if(!tree.connected(e.src, e.dst)){
      tree.link(e.src, e.dst, std::pair<int, Edge>(weight, e));
      continue;
    }

    std::pair<int, Edge> max_weight_edge = tree.path_query(e.src, e.dst);
    if(max_weight_edge.first > weight){
      Edge to_del = max_weight_edge.second; tree.cut(to_del.src, to_del.dst);
      tree.link(e.src, e.dst, std::pair<int, Edge>(weight, e));
    }
  }
}
*/
double incremental_conn_benchmark(parlay::sequence<parlay::sequence<vertex_t>>& G, long seed = -1){
  return 0.0;
}
};

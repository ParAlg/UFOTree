#include "util.h"
#include "../include/ternarized_tree.h"
#include "../include/ufo_tree.h"
#include "../include/topology_tree.h"
#include "rc_tree.h"
#include <limits>

template <typename aug_t>
void print_tree(RCTree<aug_t> *tree){
  int max_round = 0;
  for(int round : tree->round_contracted){
    max_round = std::max(max_round, round);
  }

  for(int round = 0; round <= max_round; round++){
    std::cout << "Round: " << round << "\n";
    for(int i = 0; i < tree->n; i++){
      if(tree->contraction_tree[i].size() >= round + 1 && tree->contraction_tree[i][round] != nullptr){
        std::cout << "Vertex: " << i << " :[";
        for(int j = 0; j < DEGREE_BOUND; j++){
          if(tree->contraction_tree[i][round][j] != nullptr){
            //std::cout << "HERE!";
            auto cluster = tree->contraction_tree[i][round][j];
            std::cout << "[";
            for(int k = 0; k < cluster->bv_size(); k++){
              std::cout << cluster->boundary_vertexes[k];
              if(k < cluster->bv_size() - 1){std::cout << ", ";}
            }
            std::cout << "]";
          }
        }
        std::cout << "] \n";
      }
    }
  } 
}

template<typename DynamicTree>
double incremental_MSF_benchmark(std::vector<std::pair<int, Edge>>& edges, int n){
  Edge garbage; garbage.src = MAX_VERTEX_T; garbage.dst = MAX_VERTEX_T;
  std::pair<int, Edge> id(std::numeric_limits<int>::min(), garbage);
  RCTree<std::pair<int, Edge>> tree(n, PATH, [] (std::pair<int, Edge> a, std::pair<int, Edge> b){return a.first > b.first ? a : b;}, id, id);

  for(auto edge_pair : edges){
    int weight = edge_pair.first;
    Edge e = edge_pair.second;

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
  print_tree(&tree);
  return 0;
}

int main() {
    Edge e; e.src = 0; e.dst = 1; std::pair<int,Edge> p(1, e);
    Edge e2; e2.src = 1;e2.dst = 2; std::pair<int,Edge> p2(3, e2);
    Edge e3; e3.src = 0; e3.dst = 2; std::pair<int,Edge> p3(2, e3);
    std::vector<std::pair<int, Edge>> v = {p, p2, p3};
    incremental_MSF_benchmark<RCTree<std::pair<int, Edge>>>(v, 4); 
}

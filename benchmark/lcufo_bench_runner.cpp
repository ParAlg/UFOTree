#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include <parlay/internal/get_time.h>
#include "util.h"
#include "ufo_tree.h"
#include "../baselines/dynamic_trees/link_cut_tree/include/link_cut_tree.hpp"
#include "types.h"

int N = 100000;

void run_64ary_ufo(std::vector<Update> updates);
void run_64ary_lc(std::vector<Update> updates);
std::vector<Update> k_ary_tree_gen(vertex_t n, long seed);

int main(){
    std::vector<Update> updates = k_ary_tree_gen(N,1); 
    run_64ary_ufo(updates);
    run_64ary_lc(updates);
}

void run_64ary_ufo(std::vector<Update> updates){
    UFOTree<int, int> t(N);
    for(auto update : updates){
        t.link(update.edge.src, update.edge.dst, 0);
    }
}

void run_64ary_lc(std::vector<Update> updates){
    link_cut_tree::LinkCutTreeInt t(N);
    for(auto update : updates){
        t.link(update.edge.src, update.edge.dst, 0);
    }
}

std::vector<Update> k_ary_tree_gen(vertex_t n, long seed) {
    vertex_t k = 64;
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(seed);
    parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
    ids = parlay::random_shuffle(ids, parlay::random(rand()));

    for (vertex_t i = 1; i < n; i++)
        edges.push_back({ids[i],ids[(i-1)/k]});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({INSERT,edge});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({DELETE,edge}); 
    
    return updates;
}
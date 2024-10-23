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

};
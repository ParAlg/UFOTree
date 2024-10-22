#pragma once
#include <vector>
#include <iostream>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include <parlay/internal/get_time.h>
#include "types.h"
#include "graph_utils.h"


namespace graph_benchmark {

std::vector<Update> stream_file_BFS_tree_benchmark(std::string file_name) {
    auto G = graph_utils::break_sym_graph_from_bin(file_name);
    auto E = graph_utils::BFS_forest(G);
    graph_utils::print_graph_stats(G);
    auto insertions = parlay::map(parlay::random_shuffle(E), [&](auto e) -> Update {return {INSERT, e};});
    auto deletions = parlay::map(parlay::random_shuffle(E), [&](auto e) -> Update {return {DELETE, e};});
    return parlay::append(insertions, deletions).to_vector();
}

};
#pragma once
#include <vector>
#include <iostream>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include <parlay/internal/get_time.h>
#include "types.h"


namespace ufo {

namespace parallel_dynamic_tree_benchmark {

// Converts the sequence of updates into ~n/k batches of size k
static std::vector<UpdateBatch> convert_updates_to_batches(const std::vector<Update>& updates, size_t batch_size) {
    std::vector<UpdateBatch> batches;
    if (updates.empty()) return batches;

    UpdateType current_type = updates[0].type;
    UpdateBatch current_batch;
    current_batch.type = current_type;
    for (const auto& update : updates) {
        if (update.type != current_type || current_batch.edges.size() >= batch_size) {
            batches.push_back(current_batch);
            current_batch.edges.clear();
            current_type = update.type;
            current_batch.type = current_type;
        }
        current_batch.edges.push_back({(int)update.edge.src, (int)update.edge.dst});
    }
    if (!current_batch.edges.empty()) {
        batches.push_back(current_batch);
    }
    return batches;
}

// Returns the time in seconds to perform all of the updates
template <typename DynamicTree>
double get_update_speed(vertex_t n, vertex_t k, std::vector<std::vector<UpdateBatch>>& update_sequences) {
    parlay::internal::timer my_timer("");
    for (auto updates : update_sequences) {
        DynamicTree tree(n, k);
        my_timer.start();
        for (auto batch : updates) {
            if (batch.type == INSERT) tree.batch_link(batch.edges);
            else tree.batch_cut(batch.edges);
        }
        my_timer.stop();
    }
    return my_timer.total_time()/update_sequences.size();
}

template <typename DynamicTree>
double get_update_speed_with_rand_edge_weights(vertex_t n, vertex_t k, std::vector<std::vector<UpdateBatchWithWeights>>& update_sequences)
{
    parlay::internal::timer my_timer("");
    for (auto updates : update_sequences){
        DynamicTree tree(n, k);
        my_timer.start();
        for (auto batch : updates){
            if (batch.type == INSERT)
                tree.batch_link(batch.insert_edges);
            else
                tree.batch_cut(batch.delete_edges);
        }
        my_timer.stop();
    }
    return my_timer.total_time() / update_sequences.size();
}
}

}

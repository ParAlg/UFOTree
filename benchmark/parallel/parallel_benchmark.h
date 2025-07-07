#pragma once
#include <vector>
#include <iostream>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include <parlay/internal/get_time.h>
#include "types.h"
#include <../spaa_rc_tree/RCtrees/ternarizer.h>

namespace parallel_dynamic_tree_benchmark {

// Returns the time in seconds to perform all of the updates
template <typename DynamicTree>
double get_update_speed(vertex_t n, vertex_t k, std::vector<std::vector<Update>> update_sequences) {
    parlay::internal::timer my_timer("");
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
    return my_timer.total_time()/update_sequences.size();
}

}

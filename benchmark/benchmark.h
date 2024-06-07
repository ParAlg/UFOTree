#include <vector>
#include <iostream>
#include <parlay/internal/get_time.h>
#include "../include/types.h"


namespace dynamic_tree_benchmark {

template <typename DynamicTree>
void perform_sequential_updates(DynamicTree* tree, std::vector<Update> updates) {
    parlay::internal::timer timer;
    timer.start();
    std::cout << "Doing " << updates.size() << " updates..." << std::endl;
    for (Update update : updates) {
        if (update.type == INSERT) {
            tree->link(update.edge.src, update.edge.dst);
        } else if (update.type == DELETE) {
            tree->cut(update.edge.src, update.edge.dst);
        } else {
            std::cerr << "Invalid update type: " << update.type << std::endl;
            std::abort();
        }
    }
    timer.next("Time for all updates: ");
}

template <typename DynamicTree>
void run_sequential_benchmark(vertex_t n) {
    DynamicTree tree(n);
    std::vector<Update> updates;
    for (vertex_t i = 0; i < n-1; i++)
        updates.push_back({INSERT,{i,i+1}});
    perform_sequential_updates<DynamicTree>(&tree, updates);
}

}
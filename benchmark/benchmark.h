#include <vector>
#include <iostream>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include <parlay/internal/get_time.h>
#include "../include/types.h"


namespace dynamic_tree_benchmark {

template <typename DynamicTree>
void perform_sequential_updates(DynamicTree* tree, std::vector<Update> updates) {
    parlay::internal::timer timer("");
    timer.start();
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
    timer.next("");
}

template <typename DynamicTree>
void incremental_linked_list_benchmark(vertex_t n) {
    DynamicTree tree(n);
    std::vector<Update> updates;
    for (vertex_t i = 0; i < n-1; i++)
        updates.push_back({INSERT,{i,i+1}});
    perform_sequential_updates<DynamicTree>(&tree, updates);
}

template <typename DynamicTree>
void random_degree3_benchmark(vertex_t n) {
    DynamicTree tree(n);
    std::vector<Update> updates;
    srand(time(NULL));
    auto seed = rand();
    parlay::sequence<Edge> edges;
    std::vector<int> vertex_degrees(n,0);
    while (edges.size() < n-1) {
        vertex_t u = edges.size()+1;
        vertex_t v = rand() % u;
        if (vertex_degrees[v] >= 3) continue;
        edges.push_back({u,v});
        vertex_degrees[u]++;
        vertex_degrees[v]++;
    }
    edges = parlay::random_shuffle(edges, parlay::random(seed));
    for (auto edge : edges) updates.push_back({INSERT,edge});
    edges = parlay::random_shuffle(edges, parlay::random(seed));
    for (auto edge : edges) updates.push_back({DELETE,edge});
    perform_sequential_updates<DynamicTree>(&tree, updates);
}

template <typename DynamicTree>
void random_unbounded_benchmark(vertex_t n) {
    DynamicTree tree(n);
    std::vector<Update> updates;
    srand(time(NULL));
    auto seed = rand();
    parlay::sequence<Edge> edges;
    while (edges.size() < n-1) {
        vertex_t u = edges.size()+1;
        vertex_t v = rand() % u;
        edges.push_back({u,v});
    }
    edges = parlay::random_shuffle(edges, parlay::random(seed));
    for (auto edge : edges) updates.push_back({INSERT,edge});
    edges = parlay::random_shuffle(edges, parlay::random(seed));
    for (auto edge : edges) updates.push_back({DELETE,edge});
    perform_sequential_updates<DynamicTree>(&tree, updates);
}

template <typename DynamicTree>
void preferential_attachment_benchmark(vertex_t n) {
    DynamicTree tree(n);
    std::vector<Update> updates;
    srand(time(NULL));
    auto seed = rand();
    parlay::sequence<Edge> edges;
    std::vector<int> vertex_degrees(n,0);
    while (edges.size() < n-1) {
        vertex_t u = edges.size()+1;
        vertex_t v = 0;
        int total_degree = 2*edges.size();
        if (total_degree > 0) { // preferential attachment
            int x = rand() % total_degree;
            int degree_sum = vertex_degrees[0];
            while (x >= degree_sum)
                degree_sum += vertex_degrees[++v];
        }
        edges.push_back({u,v});
        vertex_degrees[u]++;
        vertex_degrees[v]++;
    }
    edges = parlay::random_shuffle(edges, parlay::random(seed));
    for (auto edge : edges) updates.push_back({INSERT,edge});
    edges = parlay::random_shuffle(edges, parlay::random(seed));
    for (auto edge : edges) updates.push_back({DELETE,edge});
    perform_sequential_updates<DynamicTree>(&tree, updates);
}

}
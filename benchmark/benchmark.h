#include <vector>
#include <iostream>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include <parlay/internal/get_time.h>
#include "../include/types.h"


namespace dynamic_tree_benchmark {

// Returns the average time in seconds to perform all of the updates
template <typename DynamicTree>
double get_update_speed(vertex_t n, std::function<std::vector<Update>(vertex_t)> update_generator, int num_trials=1) {
    parlay::internal::timer my_timer("");
    for (int trial = 0; trial < num_trials; trial++) {
        DynamicTree tree(n);
        std::vector<Update> updates = update_generator(n);
        my_timer.start();
        for (Update update : updates) {
            if (update.type == INSERT) {
                tree.link(update.edge.src, update.edge.dst);
            } else if (update.type == DELETE) {
                tree.cut(update.edge.src, update.edge.dst);
            } else {
                std::cerr << "Invalid update type: " << update.type << std::endl;
                std::abort();
            }
        }
        my_timer.stop();
    }
    return my_timer.total_time()/num_trials;
}

// TODO: query speed
template <typename DynamicTree>
double get_query_speed(vertex_t n, std::vector<Update> updates) {
    return 0;
}

// Returns space in bytes usage just before the first deletion
template <typename DynamicTree>
size_t get_peak_space(vertex_t n, std::vector<Update> updates) {
    DynamicTree tree(n);
    size_t space;
    bool first_delete = true;
    for (Update update : updates) {
        if (update.type == INSERT) {
            tree->link(update.edge.src, update.edge.dst);
        } else if (update.type == DELETE) {
            if (first_delete) {
                space = tree->space();
                first_delete = false;
            }
            tree->cut(update.edge.src, update.edge.dst);
        } else {
            std::cerr << "Invalid update type: " << update.type << std::endl;
            std::abort();
        }
    }
    return space;
}

/* ==============================================================================
Each benchmark class returns a randomly ordered list of updates for its test case
============================================================================== */

std::vector<Update> linked_list_benchmark(vertex_t n) {
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(time(NULL));
    parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
    ids = parlay::random_shuffle(ids, parlay::random(rand()));

    for (vertex_t i = 0; i < n-1; i++)
        edges.push_back({ids[i],ids[i+1]});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({INSERT,edge});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({DELETE,edge});

    return updates;
}

std::vector<Update> binary_tree_benchmark(vertex_t n) {
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(time(NULL));
    parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
    ids = parlay::random_shuffle(ids, parlay::random(rand()));

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        edges.push_back({ids[i],ids[2*i+1]});
        edges.push_back({ids[i],ids[2*i+2]});
    }
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({INSERT,edge});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({DELETE,edge});
    
    return updates;
}

std::vector<Update> k_ary_tree_benchmark(vertex_t n) {
    vertex_t k = 64;
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(time(NULL));
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

std::vector<Update> star_benchmark(vertex_t n) {
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(time(NULL));
    parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
    ids = parlay::random_shuffle(ids, parlay::random(rand()));

    for (vertex_t i = 1; i < n; i++)
        edges.push_back({ids[0],ids[i]});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({INSERT,edge});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({DELETE,edge});
    
    return updates;
}

std::vector<Update> random_degree3_benchmark(vertex_t n) {
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(time(NULL));
    parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
    ids = parlay::random_shuffle(ids, parlay::random(rand()));

    std::vector<int> vertex_degrees(n,0);
    while (edges.size() < n-1) {
        vertex_t u = edges.size()+1;
        vertex_t v = rand() % u;
        if (vertex_degrees[v] >= 3) continue;
        edges.push_back({ids[u],ids[v]});
        vertex_degrees[u]++;
        vertex_degrees[v]++;
    }
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({INSERT,edge});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({DELETE,edge});
    
    return updates;
}

std::vector<Update> random_unbounded_benchmark(vertex_t n) {
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(time(NULL));
    parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
    ids = parlay::random_shuffle(ids, parlay::random(rand()));

    while (edges.size() < n-1) {
        vertex_t u = edges.size()+1;
        vertex_t v = rand() % u;
        edges.push_back({ids[u],ids[v]});
    }
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({INSERT,edge});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({DELETE,edge});
    
    return updates;
}

std::vector<Update> preferential_attachment_benchmark(vertex_t n) {
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(time(NULL));
    parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
    ids = parlay::random_shuffle(ids, parlay::random(rand()));

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
        edges.push_back({ids[u],ids[v]});
        vertex_degrees[u]++;
        vertex_degrees[v]++;
    }
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({INSERT,edge});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({DELETE,edge});

    return updates;
}

}

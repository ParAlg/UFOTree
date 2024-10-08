#pragma once
#include <vector>
#include <iostream>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include <parlay/internal/get_time.h>
#include "types.h"


namespace dynamic_tree_benchmark {

// Returns the time in seconds to perform all of the updates
template <typename DynamicTree>
double get_update_speed(vertex_t n, std::vector<std::vector<Update>> update_sequences) {
    parlay::internal::timer my_timer("");
    for (auto updates : update_sequences) {
        DynamicTree tree(n);
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
    return my_timer.total_time()/update_sequences.size();
}

// TODO: query speed
template <typename DynamicTree>
double get_query_speed(vertex_t n, std::vector<std::vector<Update>> update_sequences, std::vector<std::vector<Query>> query_sequences) {
    parlay::internal::timer my_timer("");
    for (int i = 0; i < update_sequences.size(); ++i) {
        auto updates = update_sequences[i];
        auto queries = query_sequences[i];
        DynamicTree tree(n);
        bool first_delete = true;
        for (Update update : updates) {
            if (update.type == INSERT) {
                tree.link(update.edge.src, update.edge.dst);
            } else if (update.type == DELETE) {
                if (first_delete) {
                    my_timer.start();
                    // for (auto query : queries) tree.path_query(query.u, query.v);
                    my_timer.stop();
                    first_delete = false;
                }
                tree.cut(update.edge.src, update.edge.dst);
            } else {
                std::cerr << "Invalid update type: " << update.type << std::endl;
                std::abort();
            }
        }
    }
    return my_timer.total_time()/query_sequences.size();
}

// Returns space in bytes used just before the first deletion
template <typename DynamicTree>
size_t get_peak_space(vertex_t n, std::vector<std::vector<Update>> update_sequences) {
    size_t total_space = 0;
    for (auto updates : update_sequences) {
        DynamicTree tree(n);
        bool first_delete = true;
        for (Update update : updates) {
            if (update.type == INSERT) {
                tree.link(update.edge.src, update.edge.dst);
            } else if (update.type == DELETE) {
                if (first_delete) {
                    total_space += tree.space();
                    first_delete = false;
                }
                tree.cut(update.edge.src, update.edge.dst);
            } else {
                std::cerr << "Invalid update type: " << update.type << std::endl;
                std::abort();
            }
        }
    }
    return total_space/update_sequences.size();
}

/* ==============================================================================
Each benchmark class returns a randomly ordered list of updates for its test case
============================================================================== */

std::vector<Update> linked_list_benchmark(vertex_t n, long seed) {
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(seed);
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

std::vector<Update> binary_tree_benchmark(vertex_t n, long seed) {
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(seed);
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

std::vector<Update> k_ary_tree_benchmark(vertex_t n, long seed) {
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

std::vector<Update> star_benchmark(vertex_t n, long seed) {
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(seed);
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

std::vector<Update> random_degree3_benchmark(vertex_t n, long seed) {
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(seed);
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

std::vector<Update> random_unbounded_benchmark(vertex_t n, long seed) {
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(seed);
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

std::vector<Update> preferential_attachment_benchmark(vertex_t n, long seed) {
    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(seed);
    parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
    ids = parlay::random_shuffle(ids, parlay::random(rand()));

    while (edges.size() < n-1) {
        vertex_t u = edges.size()+1;
        vertex_t v = 0;
        if (edges.size() > 0) {
            auto rand_edge = edges[rand() % edges.size()];
            v = rand()%2 ? rand_edge.src : rand_edge.dst;
        }
        edges.push_back({u,v});
    }
    for (int i = 0; i < edges.size(); i++) {
        auto edge = edges[i];
        edges[i] = {ids[edge.src],ids[edge.dst]};
    }
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({INSERT,edge});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({DELETE,edge});

    return updates;
}

std::vector<Query> random_query_generator(vertex_t n, vertex_t num_queries) {
    std::vector<Query> queries;
    srand(time(NULL));
    for (int i = 0; i < num_queries; ++i)
        queries.push_back({rand()%n, rand()%n});
    return queries;
}

}

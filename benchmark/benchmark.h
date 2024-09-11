#include <vector>
#include <iostream>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include <parlay/internal/get_time.h>
#include "../include/types.h"

vertex_t n_list[] = {10,100};

namespace dynamic_tree_benchmark {

#define COLLECT_SPACE


static parlay::internal::timer timer("");

template <typename DynamicTree>
void perform_sequential_updates(DynamicTree* tree, std::vector<Update> updates, size_t* space, bool init_calc=false) {
    #ifdef COLLECT_SPACE
      if(!init_calc){
        size_t initial_space = tree->space();
        std::cout << "\nInit Space: " << initial_space << " Bytes" << std::endl;
      }
    #endif
    timer.start();
    bool first_delete = true;
    for (Update update : updates) {
        if (update.type == INSERT) {
            tree->link(update.edge.src, update.edge.dst);
        } else if (update.type == DELETE) {
            #ifdef COLLECT_SPACE
            if (first_delete) {
                timer.stop();
                *space = tree->space();
                first_delete = false;
                timer.start();
            } 
            #endif
            tree->cut(update.edge.src, update.edge.dst);
        } else {
            std::cerr << "Invalid update type: " << update.type << std::endl;
            std::abort();
        }
    }
}

template <typename DynamicTree>
void linked_list_benchmark(vertex_t n) {
    DynamicTree tree(n);
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

    size_t space_used = 0;
    timer.reset();
    timer.start();
    perform_sequential_updates<DynamicTree>(&tree, updates, &space_used);
    std::cout << "Total Time";
    timer.next("");
    #ifdef COLLECT_SPACE 
      std::cout << "Peak Space: " << space_used << " Bytes" << std::endl;
    #endif
}

template <typename DynamicTree>
void binary_tree_benchmark(vertex_t n) {
    DynamicTree tree(n);
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
    
    size_t space_used = 0;
    timer.reset();
    perform_sequential_updates<DynamicTree>(&tree, updates, &space_used);
    std::cout << "Total Time";
    timer.next("");
    #ifdef COLLECT_SPACE 
      std::cout << "Peak Space: " << space_used << " Bytes" << std::endl;
    #endif
}

template <typename DynamicTree>
void k_ary_tree_benchmark(vertex_t n, vertex_t k = 64) {
    DynamicTree tree(n);
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
    
    size_t space_used = 0;
    timer.reset();
    perform_sequential_updates<DynamicTree>(&tree, updates, &space_used);
    std::cout << "Total Time";
    timer.next("");
     #ifdef COLLECT_SPACE 
      std::cout << "Peak Space: " << space_used << " Bytes" << std::endl;
    #endif
}

template <typename DynamicTree>
void star_benchmark(vertex_t n) {
    DynamicTree tree(n);
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
    size_t space_used = 0;
    timer.reset();
    perform_sequential_updates<DynamicTree>(&tree, updates, &space_used);
    std::cout << "Total Time";
    timer.next("");
    #ifdef COLLECT_SPACE 
      std::cout << "Peak Space: " << space_used << " Bytes" << std::endl;
    #endif
}

template <typename DynamicTree>
void random_degree3_benchmark(vertex_t n) {
    int num_trials = 5;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    long long space_total = 0;
    size_t space_used = 0;
    for (int trial = 0; trial < num_trials; trial++) {
        DynamicTree tree(n);
        std::vector<Update> updates;
        parlay::sequence<Edge> edges;
        std::vector<int> vertex_degrees(n,0);
        auto seed = seeds[trial];
        srand(seed);
        parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
        ids = parlay::random_shuffle(ids, parlay::random(rand()));

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
        
        timer.reset();
        if(trial == 0){
          perform_sequential_updates<DynamicTree>(&tree, updates, &space_used);
        } else {
          perform_sequential_updates(&tree, updates, &space_used, true);
        }
        timer.stop();
        space_total += space_used;
    }
    std::cout << "Total Time: " << timer.total_time()/num_trials << std::endl;
    #ifdef COLLECT_SPACE
      std::cout << "Peak Space: " << space_total/num_trials << " Bytes" << std::endl;
    #endif
    
}

template <typename DynamicTree>
void random_unbounded_benchmark(vertex_t n) {
    int num_trials = 5;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    parlay::internal::timer timer("");
    long long space_total = 0;

    for (int trial = 0; trial < num_trials; trial++) {
        DynamicTree tree(n);
        std::vector<Update> updates;
        parlay::sequence<Edge> edges;
        auto seed = seeds[trial];
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
        
        size_t space_used = 0;
        timer.reset();
        if(trial == 0){
          perform_sequential_updates<DynamicTree>(&tree, updates, &space_used);
        } else {
          perform_sequential_updates(&tree, updates, &space_used, true);
        }
        timer.stop();
        space_total += space_used;
    }
    std::cout << "Total Time: " << timer.total_time()/num_trials << std::endl; 
    #ifdef COLLECT_SPACE
      std::cout << "Peak Space: " << space_total/num_trials << " Bytes" << std::endl;
    #endif
}

template <typename DynamicTree>
void preferential_attachment_benchmark(vertex_t n) {
    int num_trials = 5;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    parlay::internal::timer timer("");
    long long space_total =  0;

    for (int trial = 0; trial < num_trials; trial++) {
        DynamicTree tree(n);
        std::vector<Update> updates;
        parlay::sequence<Edge> edges;
        std::vector<int> vertex_degrees(n,0);
        auto seed = seeds[trial];
        srand(seed);
        parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
        ids = parlay::random_shuffle(ids, parlay::random(rand()));

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
        size_t space_used = 0; 
        timer.reset();
        if(trial == 0){
          perform_sequential_updates<DynamicTree>(&tree, updates, &space_used);
        } else {
          perform_sequential_updates(&tree, updates, &space_used, true);
        }
        timer.stop();
        space_total += space_used;
    }
    std::cout << "Total Time: " << timer.total_time()/num_trials << std::endl;
    #ifdef COLLECT_SPACE
      std::cout << "Peak Space: " << space_total/num_trials << " Bytes" << std::endl;
    #endif
}

}

#include <gtest/gtest.h>
#include <unordered_set>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include "../include/parallel_topology_tree.h"


template<typename aug_t>
bool ParallelTopologyTree<aug_t>::is_valid() {
    std::unordered_set<ParallelTopologyCluster<aug_t>*> clusters;
    std::unordered_set<ParallelTopologyCluster<aug_t>*> next_clusters;
    for (auto leaf : leaves) // Ensure that every pair of incident vertices are in the same component
        for (auto neighbor : leaf->neighbors) // This ensures all connectivity is correct by transitivity
            if (neighbor && leaf->get_root() != neighbor->get_root()) return false;
    for (int i = 0; i < this->leaves.size(); i++) clusters.insert(leaves[i]);
    while (!clusters.empty()) {
        for (auto cluster : clusters) {
            for (auto neighbor : cluster->neighbors) // Ensure all neighbors also point back
                if (neighbor && !neighbor->contains_neighbor(cluster)) return false;
            if (!cluster->contracts()) { // Ensure maximality of contraction
                if (cluster->get_degree() == 1) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && !neighbor->contracts()) return false;
                } else if (cluster->get_degree() == 2) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && !neighbor->contracts() && neighbor->get_degree() < 3) return false;
                } else if (cluster->get_degree() == 3) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && !neighbor->contracts() && neighbor->get_degree() < 2) return false;
                }
            }
            if (cluster->parent) next_clusters.insert(cluster->parent); // Get next level
        }
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
    return true;
}

template<typename aug_t>
int ParallelTopologyTree<aug_t>::get_height(vertex_t v) {
    int height = 0;
    ParallelTopologyCluster<aug_t>* curr = &leaves[v];
    while (curr) {
        height++;
        curr = curr->parent;
    }
    return height;
}

template<typename aug_t>
void ParallelTopologyTree<aug_t>::print_tree() {
    std::multimap<ParallelTopologyCluster<aug_t>*, ParallelTopologyCluster<aug_t>*> clusters;
    std::multimap<ParallelTopologyCluster<aug_t>*, ParallelTopologyCluster<aug_t>*> next_clusters;
    std::cout << "========================= LEAVES =========================" << std::endl;
    std::unordered_map<ParallelTopologyCluster<aug_t>*, vertex_t> vertex_map;
    for (int i = 0; i < this->leaves.size(); i++) vertex_map.insert({leaves[i], i});
    for (int i = 0; i < this->leaves.size(); i++) clusters.insert({leaves[i]->parent, leaves[i]});
    for (auto entry : clusters) {
        auto leaf = entry.second;
        auto parent = entry.first;
        std::cout << "VERTEX " << vertex_map[leaf] << "\t " << leaf << " Parent " << parent << " Neighbors: ";
        for (auto neighbor : leaf->neighbors) if (neighbor) std::cout << vertex_map[neighbor] << " ";
        std::cout << std::endl;
        bool in_map = false;
        for (auto entry : next_clusters) if (entry.second == parent) in_map = true;
        if (parent && !in_map) next_clusters.insert({parent->parent, parent});
    }
    clusters.swap(next_clusters);
    next_clusters.clear();
    while (!clusters.empty()) {
        std::cout << "======================= NEXT LEVEL =======================" << std::endl;
        for (auto entry : clusters) {
            auto cluster = entry.second;
            auto parent = entry.first;
            std::cout << "Cluster: " << cluster << " Parent: " << parent << std::endl;
            bool in_map = false;
            for (auto entry : next_clusters) if (entry.second == parent) in_map = true;
            if (parent && !in_map) next_clusters.insert({parent->parent, parent});
        }
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
}

TEST(ParallelTopologyTreeSuite, batch_incremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    vertex_t k = 1;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    ParallelTopologyTree<int> tree(n, qt, f, 0, 0);

    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(time(NULL));
    parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
    ids = parlay::random_shuffle(ids, parlay::random(rand()));
    for (vertex_t i = 0; i < n-1; i++)
        edges.push_back({ids[i],ids[i+1]});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({INSERT,edge});

    Edge batch[k];
    vertex_t len = 0;
    for (auto update : updates) {
        batch[len++] = update.edge;
        if (len == k) {
            tree.batch_link(batch, len);
            len = 0;
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
        }
    }
    tree.batch_link(batch, len);
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
}

TEST(ParallelTopologyTreeSuite, decremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 1000000;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int{return x + y;};
        ParallelTopologyTree<int> tree(n, qt, f, 0, 0);
        std::pair<vertex_t, vertex_t> edges[n-1];

        auto seed = seeds[trial];
        srand(seed);
        int links = 0;
        std::vector<int> vertex_degrees(n,0);
        while (links < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (vertex_degrees[u] >= 3) continue;
            if (vertex_degrees[v] >= 3) continue;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                edges[links++] = {u,v};
                vertex_degrees[u]++;
                vertex_degrees[v]++;
            }
        }
        for (auto edge : edges) {
            auto u = edge.first;
            auto v = edge.second;
            tree.cut(u,v);
            ASSERT_FALSE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " connected.";
            // ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << u << " and " << v << ".";
        }
    }
}

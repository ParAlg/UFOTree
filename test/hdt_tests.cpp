#include <gtest/gtest.h>
#include <unordered_set>
#include "../include/hdt_ufo_tree.h"


bool HDTUFOTree::is_valid() {
    std::unordered_set<HDTUFOCluster*> clusters;
    std::unordered_set<HDTUFOCluster*> next_clusters;
    for (auto leaf : leaves) { // Ensure that every pair of incident vertices are in the same component
        for (auto neighbor : leaf.neighbors) // This ensures all connectivity is correct by transitivity
            if (neighbor && leaf.get_root() != neighbor->get_root()) return false;
        if (leaf.neighbors_set)
        for (auto neighbor : *leaf.neighbors_set) {
            if (leaf.get_root() != neighbor->get_root()) return false;
        }
    }
    for (int i = 0; i < this->leaves.size(); i++) clusters.insert(&this->leaves[i]);
    while (!clusters.empty()) {
        for (auto cluster : clusters) {
            for (auto neighbor : cluster->neighbors) // Ensure all neighbors also point back
                if (neighbor && !neighbor->contains_neighbor(cluster)) return false;
            if (cluster->neighbors_set)
            for (auto neighbor : *cluster->neighbors_set) {
                if (!neighbor->contains_neighbor(cluster)) return false;
            }
            if (cluster->get_degree() <= 3 && !cluster->contracts()) { // Ensure maximality of contraction
                if (cluster->get_degree() == 1) {
                    if (cluster->get_neighbor()->get_degree() > 2) return false;
                    else if (!cluster->get_neighbor()->contracts()) return false;
                } else if (cluster->get_degree() == 2) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && neighbor->get_degree() < 3 && !neighbor->contracts()) return false;
                } else if (cluster->get_degree() >= 3) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && neighbor->get_degree() < 2) return false;
                    if (cluster->neighbors_set)
                    for (auto neighbor : *cluster->neighbors_set)
                        if (neighbor->get_degree() < 2) return false;
                }
            }
            if (cluster->parent) next_clusters.insert(cluster->parent); // Get next level
        }
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
    return true;
}

void HDTUFOTree::print_tree() {
    std::multimap<HDTUFOCluster*, HDTUFOCluster*> clusters;
    std::multimap<HDTUFOCluster*, HDTUFOCluster*> next_clusters;
    std::cout << "========================= LEAVES =========================" << std::endl;
    std::unordered_map<HDTUFOCluster*, vertex_t> vertex_map;
    for (int i = 0; i < this->leaves.size(); i++) vertex_map.insert({&leaves[i], i});
    for (int i = 0; i < this->leaves.size(); i++) clusters.insert({leaves[i].parent, &leaves[i]});
    for (auto entry : clusters) {
        auto leaf = entry.second;
        auto parent = entry.first;
        std::cout << "VERTEX " << vertex_map[leaf] << "\t " << leaf << " Parent " << parent << " Neighbors: ";
        for (auto neighbor : leaf->neighbors) if (neighbor) std::cout << vertex_map[neighbor] << " ";
        if (leaf->neighbors_set)
        for (auto neighbor : *leaf->neighbors_set) std::cout << vertex_map[neighbor] << " ";
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

TEST(HDTUFOTreeSuite, incremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    HDTUFOTree tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.AddEdge({i,i+1});
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " and " << i+1 << ".";
    }
}

TEST(HDTUFOTreeSuite, incremental_binarytree_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    HDTUFOTree tree(n);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.AddEdge({i,2*i+1});
        tree.AddEdge({i,2*i+2});
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " with children.";
    }
    if (n%2 == 0) tree.AddEdge({(n-1)/2,n-1});
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all links.";
}

TEST(HDTUFOTreeSuite, incremental_star_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    HDTUFOTree tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.AddEdge({0,i+1});
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << 0 << " and " << i+1 << ".";
    }
}

TEST(HDTUFOTreeSuite, incremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int{return x + y;};
        HDTUFOTree tree(n);

        auto seed = seeds[trial];
        srand(seed);
        int links = 0;
        while (links < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (u != v && !tree.IsConnected(u,v)) {
                tree.AddEdge({u,v});
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << u << " and " << v << ".";
                links++;
            }
        }
    }
}

TEST(HDTUFOTreeSuite, decremental_linkedlist_correctness_test) {
    vertex_t n = 128;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    HDTUFOTree tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.AddEdge({i,i+1});
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.DeleteEdge({i,i+1});
        ASSERT_FALSE(tree.IsConnected(i,i+1)) << "Vertex " << i << " and " << i+1 << " IsConnected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " and " << i+1 << ".";
    }
}

TEST(HDTUFOTreeSuite, decremental_binarytree_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    HDTUFOTree tree(n);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.AddEdge({i,2*i+1});
        tree.AddEdge({i,2*i+2});
    }
    if (n%2 == 0) tree.AddEdge({(n-1)/2,n-1});
    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.DeleteEdge({i,2*i+1});
        tree.DeleteEdge({i,2*i+2});
        ASSERT_FALSE(tree.IsConnected(i,2*i+1)) << "Vertex " << i << " and " << 2*i+1 << " IsConnected.";
        ASSERT_FALSE(tree.IsConnected(i,2*i+2)) << "Vertex " << i << " and " << 2*i+2 << " IsConnected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " from children.";
    }
    if (n%2 == 0) tree.DeleteEdge({(n-1)/2,n-1});
    ASSERT_FALSE(tree.IsConnected((n-1)/2,n-1)) << "Vertex " << (n-1)/2 << " and " << n-1 << " IsConnected.";
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all cuts.";
}

TEST(HDTUFOTreeSuite, decremental_star_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    HDTUFOTree tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.AddEdge({0,i+1});
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " and " << i+1 << ".";
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.DeleteEdge({0,i+1});
        ASSERT_FALSE(tree.IsConnected(0,i+1)) << "Vertex " << 0 << " and " << i+1 << " IsConnected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << 0 << " and " << i+1 << ".";
    }
}

TEST(HDTUFOTreeSuite, decremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int{return x + y;};
        HDTUFOTree tree(n);
        std::pair<vertex_t, vertex_t> edges[n-1];

        auto seed = seeds[trial];
        srand(seed);
        int links = 0;
        while (links < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (u != v && !tree.IsConnected(u,v)) {
                tree.AddEdge({u,v});
                edges[links++] = {u,v};
            }
        }
        for (auto edge : edges) {
            auto u = edge.first;
            auto v = edge.second;
            tree.DeleteEdge({u,v});
            ASSERT_FALSE(tree.IsConnected(u,v)) << "Vertex " << u << " and " << v << " IsConnected.";
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << u << " and " << v << ".";
        }
    }
}

TEST(HDTUFOTreeSuite, size_augmentation_test) {
    vertex_t n = 256;
    HDTUFOTree tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.AddEdge({i,i+1});
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.DeleteEdge({i,i+1});
        ASSERT_EQ(tree.GetSizeOfTree(i), i+1);
        ASSERT_EQ(tree.GetSizeOfTree(i+1), n-(i+1));
        tree.AddEdge({i,i+1});
    }
}

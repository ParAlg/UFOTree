#include <gtest/gtest.h>
#include "../include/ufo_tree.h"


template<typename aug_t>
bool UFOTree<aug_t>::is_valid() {
    std::unordered_set<UFOCluster<aug_t>*> clusters;
    std::unordered_set<UFOCluster<aug_t>*> next_clusters;
    for (auto leaf : leaves) // Ensure that every pair of incident vertices are in the same component
        for (auto entry : leaf.neighbors) // This ensures all connectivity is correct by transitivity
            if (leaf.get_root() != entry.first->get_root()) return false;
    for (int i = 0; i < this->leaves.size(); i++) clusters.insert(&this->leaves[i]);
    while (!clusters.empty()) {
        for (auto cluster : clusters) {
            for (auto entry : cluster->neighbors) // Ensure all neighbors also point back
                if (!entry.first->contains_neighbor(cluster)) return false;
            if (!cluster->contracts()) { // Ensure maximality of contraction
                if (cluster->get_degree() == 1) return false;
                else if (cluster->get_degree() == 2) {
                    for (auto entry : cluster->neighbors)
                        if (!entry.first->contracts() && entry.first->get_degree() < 3) return false;
                } else if (cluster->get_degree() >= 3) {
                    for (auto entry : cluster->neighbors)
                        if (entry.first->get_degree() < 2) return false;
                }
            }
            if (cluster->parent) next_clusters.insert(cluster->parent); // Get next level
        }
        // Maximality should ensure this, but we leave the test as a sanity check
        if (6*next_clusters.size() > 5*clusters.size()) return false;
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
    return true;
}

TEST(UFOTreeSuite, incremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    UFOTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        for (vertex_t u = 0; u < i+1; u++) for (vertex_t v = u+1; v <= i+1; v++)
            ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " and " << i+1 << ".";
    }
}

TEST(UFOTreeSuite, incremental_binarytree_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    UFOTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
        for (vertex_t u = 0; u < 2*i+2; u++) for (vertex_t v = u+1; v <= 2*i+2; v++)
            ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " with children.";
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    for (vertex_t u = 0; u < n-1; u++) for (vertex_t v = u+1; v < n; v++)
        ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all links.";
}

TEST(UFOTreeSuite, incremental_star_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    UFOTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(0,i+1);
        for (vertex_t u = 0; u < i+1; u++) for (vertex_t v = u+1; v <= i+1; v++)
            ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " and " << i+1 << ".";
    }
}

TEST(UFOTreeSuite, incremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int{return x + y;};
        UFOTree<int> tree(n, qt, f, 0, 0);

        auto seed = seeds[trial];
        srand(seed);
        int links = 0;
        while (links < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << u << " and " << v << ".";
                links++;
            }
        }
    }
}

TEST(UFOTreeSuite, decremental_linkedlist_correctness_test) {
    vertex_t n = 128;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    UFOTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        std::cout << "CUTTING " << i << " " << i+1 << std::endl;
        tree.cut(i,i+1);
        for (vertex_t u = 0; u < i+1; u++) for (vertex_t v = u+1; v < n; v++)
            ASSERT_FALSE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " connected.";
        for (vertex_t u = i+1; u < n-1; u++) for (vertex_t v = u+1; v < n; v++)
            ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " and " << i+1 << ".";
    }
}

TEST(UFOTreeSuite, decremental_binarytree_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    UFOTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.cut(i,2*i+1);
        tree.cut(i,2*i+2);
        for (vertex_t u = 0; u < i+1; u++) for (vertex_t v = u+1; v < n; v++)
            ASSERT_FALSE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " connected.";
        ASSERT_FALSE(tree.connected(2*i+1,2*i+2)) << "Vertex " << 2*i+1 << " and " << 2*i+2 << " connected.";
        if (4*i+6 < n) {
            ASSERT_TRUE(tree.connected(4*i+3,4*i+4)) << "Vertex " << 4*i+3 << " and " << 4*i+4 << " not connected.";
            ASSERT_TRUE(tree.connected(4*i+5,4*i+6)) << "Vertex " << 4*i+3 << " and " << 4*i+4 << " not connected.";
        }
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " from children.";
    }
    if (n%2 == 0) tree.cut((n-1)/2,n-1);
    for (vertex_t u = 0; u < n-1; u++) for (vertex_t v = u+1; v < n; v++)
        ASSERT_FALSE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " connected.";
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all cuts.";
}

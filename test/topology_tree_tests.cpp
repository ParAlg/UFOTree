#include <gtest/gtest.h>
#include "../include/topology_tree.h"


template<typename aug_t>
bool TopologyTree<aug_t>::is_valid() {
    std::unordered_set<TopologyCluster<aug_t>*> clusters;
    std::unordered_set<TopologyCluster<aug_t>*> next_clusters;
    for (auto leaf : leaves) // Ensure that every pair of incident vertices are in the same component
        for (auto neighbor : leaf.neighbors) // This ensures all connectivity is correct by transitivity
            if (neighbor && leaf.get_root() != neighbor->get_root()) return false;
    for (int i = 0; i < this->leaves.size(); i++) clusters.insert(&this->leaves[i]);
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

TEST(TopologyTreeSuite, incremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " and " << i+1 << ".";
    }
}

TEST(TopologyTreeSuite, incremental_binarytree_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " with children.";
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all links.";
}

TEST(TopologyTreeSuite, incremental_random_correctness_test) {
    int num_trials = 10;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int{return x + y;};
        TopologyTree<int> tree(n, qt, f, 0, 0);

        auto seed = seeds[trial];
        srand(seed);
        int links = 0;
        while (links < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (tree.leaves[u].get_degree() >= 3) continue;
            if (tree.leaves[v].get_degree() >= 3) continue;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << u << " and " << v << ".";
                links++;
            }
        }
    }
}

TEST(TopologyTreeSuite, decremental_linkedlist_correctness_test) {
    vertex_t n = 128;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(i,i+1);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " and " << i+1 << ".";
    }
}

TEST(TopologyTreeSuite, decremental_binarytree_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.cut(i,2*i+1);
        tree.cut(i,2*i+2);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " from children.";
    }
    if (n%2 == 0) tree.cut((n-1)/2,n-1);
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all cuts.";
}

TEST(TopologyTreeSuite, decremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int{return x + y;};
        TopologyTree<int> tree(n, qt, f, 0, 0);
        std::pair<vertex_t, vertex_t> edges[n-1];

        auto seed = seeds[trial];
        srand(seed);
        int links = 0;
        while (links < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (tree.leaves[u].get_degree() >= 3) continue;
            if (tree.leaves[v].get_degree() >= 3) continue;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                edges[links++] = {u,v};
            }
        }
        for (auto edge : edges) {
            auto u = edge.first;
            auto v = edge.second;
            tree.cut(u,v);
            ASSERT_FALSE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " connected.";
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << u << " and " << v << ".";
        }
    }
}

TEST(TopologyTreeSuite, subtree_query_test) {
    vertex_t n = 256;
    QueryType qt = SUBTREE;
    auto f = [](int x, int y)->int{return x + y;};
    int id = 0;
    int d = 1;
    TopologyTree<int> tree(n, qt, f, id, d);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        for (vertex_t v = 0; v <= i+1; v++) {
            ASSERT_EQ(tree.subtree_query(v, (v > 0 ? v-1 : MAX_VERTEX_T)), i+2-v);
            ASSERT_EQ(tree.subtree_query(v, (v <= i ? v+1 : MAX_VERTEX_T)), v+1);
        }
    }
    for (vertex_t i = n-1; i > 0; i--) {
        tree.cut(i-1,i);
        for (vertex_t v = 0; v < i; v++) {
            ASSERT_EQ(tree.subtree_query(v, (v > 0 ? v-1 : MAX_VERTEX_T)), i-v);
            ASSERT_EQ(tree.subtree_query(v, (v < i-1 ? v+1 : MAX_VERTEX_T)), v+1);
        }
    }
}

TEST(TopologyTreeSuite, path_query_test) {
    vertex_t n = 64;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    int id = 0;
    int d = 0;
    TopologyTree<int> tree(n, qt, f, id, d);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        for (vertex_t u = 0; u < i+1; u++) for (vertex_t v = u+1; v <= i+1; v++)
            ASSERT_EQ(tree.path_query(u,v), v-u);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(i,i+1);
        for (vertex_t u = i+1; u < n-1; u++) for (vertex_t v = u+1; v < n; v++)
            ASSERT_EQ(tree.path_query(u,v), v-u);
    }
}

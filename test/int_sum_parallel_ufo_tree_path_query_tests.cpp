#include <gtest/gtest.h>
#include <parlay/sequence.h>
#include <algorithm>
#include <cstdint>
#include <queue>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>
#include "int_sum_parallel_ufo_tree.h"

using namespace ufo;

namespace {

std::pair<int, int> norm_edge(int u, int v) {
    if (u > v) std::swap(u, v);
    return {u, v};
}

bool ref_connected(const std::vector<std::unordered_set<int>>& adj, int u, int v) {
    if (u == v) return true;
    std::vector<bool> vis(adj.size(), false);
    std::queue<int> q;
    q.push(u);
    vis[u] = true;
    while (!q.empty()) {
        int x = q.front();
        q.pop();
        for (int y : adj[x]) {
            if (vis[y]) continue;
            if (y == v) return true;
            vis[y] = true;
            q.push(y);
        }
    }
    return false;
}

uint32_t ref_path_sum(const std::vector<std::unordered_set<int>>& adj, const std::vector<uint32_t>& w, int u, int v) {
    std::vector<int> parent(adj.size(), -1);
    std::queue<int> q;
    q.push(u);
    parent[u] = u;
    while (!q.empty() && parent[v] == -1) {
        int x = q.front();
        q.pop();
        for (int y : adj[x]) {
            if (parent[y] != -1) continue;
            parent[y] = x;
            q.push(y);
        }
    }
    EXPECT_NE(parent[v], -1);
    uint32_t sum = 0;
    int cur = v;
    while (cur != parent[cur]) {
        sum += w[cur];
        cur = parent[cur];
    }
    sum += w[cur];
    return sum;
}

std::vector<int> ref_path_vertices(const std::vector<std::unordered_set<int>>& adj, int u, int v) {
    std::vector<int> parent(adj.size(), -1);
    std::queue<int> q;
    q.push(u);
    parent[u] = u;
    while (!q.empty() && parent[v] == -1) {
        int x = q.front();
        q.pop();
        for (int y : adj[x]) {
            if (parent[y] != -1) continue;
            parent[y] = x;
            q.push(y);
        }
    }
    std::vector<int> path;
    if (parent[v] == -1) return path;
    int cur = v;
    while (cur != parent[cur]) {
        path.push_back(cur);
        cur = parent[cur];
    }
    path.push_back(cur);
    std::reverse(path.begin(), path.end());
    return path;
}

}  // namespace

TEST(IntSumParallelUFOTreePathQuerySuite, single_edge_query) {
    ParallelUFOTree<> tree(2, 1);
    parlay::sequence<std::pair<int, int>> links = {{0, 1}};
    tree.BatchLink(links);

    tree.UpdateWeight(0, 7);
    tree.UpdateWeight(1, 11);

    EXPECT_TRUE(tree.connected(0, 1));
    EXPECT_EQ(tree.PathQuery(0, 1), 18u);
}

TEST(IntSumParallelUFOTreePathQuerySuite, linked_list_query) {
    ParallelUFOTree<> tree(4, 2);
    parlay::sequence<std::pair<int, int>> links = {{0, 1}, {1, 2}, {2, 3}};
    tree.BatchLink(links);

    tree.UpdateWeight(0, 1);
    tree.UpdateWeight(1, 2);
    tree.UpdateWeight(2, 3);
    tree.UpdateWeight(3, 4);

    EXPECT_TRUE(tree.connected(0, 3));
    EXPECT_EQ(tree.PathQuery(0, 3), 10u);
    EXPECT_EQ(tree.PathQuery(1, 2), 5u);
}

TEST(IntSumParallelUFOTreePathQuerySuite, same_vertex_query) {
    ParallelUFOTree<> tree(3, 2);
    parlay::sequence<std::pair<int, int>> links = {{0, 1}, {1, 2}};
    tree.BatchLink(links);

    tree.UpdateWeight(0, 5);
    tree.UpdateWeight(1, 13);
    tree.UpdateWeight(2, 21);

    EXPECT_EQ(tree.PathQuery(1, 1), 13u);
}

TEST(IntSumParallelUFOTreePathQuerySuite, update_weight_changes_path_result) {
    ParallelUFOTree<> tree(3, 2);
    parlay::sequence<std::pair<int, int>> links = {{0, 1}, {1, 2}};
    tree.BatchLink(links);

    tree.UpdateWeight(0, 10);
    tree.UpdateWeight(1, 20);
    tree.UpdateWeight(2, 30);
    EXPECT_EQ(tree.PathQuery(0, 2), 60u);

    tree.UpdateWeight(1, 5);
    EXPECT_EQ(tree.PathQuery(0, 2), 45u);
}

TEST(IntSumParallelUFOTreePathQuerySuite, balanced_binary_tree_multiple_queries) {
    ParallelUFOTree<> tree(7, 4);
    parlay::sequence<std::pair<int, int>> links = {
        {0, 1}, {0, 2}, {1, 3}, {1, 4}, {2, 5}, {2, 6}
    };
    tree.BatchLink(links);

    std::vector<uint32_t> w = {5, 3, 7, 11, 13, 17, 19};
    for (int i = 0; i < (int) w.size(); ++i) tree.UpdateWeight(i, w[i]);

    EXPECT_EQ(tree.PathQuery(3, 4), 11u + 3u + 13u);
    EXPECT_EQ(tree.PathQuery(3, 6), 11u + 3u + 5u + 7u + 19u);
    EXPECT_EQ(tree.PathQuery(5, 6), 17u + 7u + 19u);
    EXPECT_EQ(tree.PathQuery(0, 6), 5u + 7u + 19u);
}

TEST(IntSumParallelUFOTreePathQuerySuite, star_tree_queries_and_updates) {
    ParallelUFOTree<> tree(8, 4);
    parlay::sequence<std::pair<int, int>> links = {
        {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}, {0, 7}
    };
    tree.BatchLink(links);

    std::vector<uint32_t> w = {100, 1, 2, 3, 4, 5, 6, 7};
    for (int i = 0; i < (int) w.size(); ++i) tree.UpdateWeight(i, w[i]);

    EXPECT_EQ(tree.PathQuery(1, 7), 1u + 100u + 7u);
    EXPECT_EQ(tree.PathQuery(2, 3), 2u + 100u + 3u);

    tree.UpdateWeight(0, 42);
    EXPECT_EQ(tree.PathQuery(1, 7), 1u + 42u + 7u);
    EXPECT_EQ(tree.PathQuery(4, 6), 4u + 42u + 6u);
}

TEST(IntSumParallelUFOTreePathQuerySuite, random_dynamic_forest_reference_check) {
    constexpr int n = 24;
    constexpr int steps = 220;
    ParallelUFOTree<> tree(n, 8);
    std::vector<std::unordered_set<int>> adj(n);
    std::unordered_set<long long> edge_set;
    std::vector<uint32_t> w(n, 1);

    std::mt19937 rng(1234567);
    std::uniform_int_distribution<int> vd(0, n - 1);
    std::uniform_int_distribution<int> wd(1, 1000);
    std::uniform_real_distribution<double> op(0.0, 1.0);

    auto edge_key = [&](int u, int v) {
        auto [a, b] = norm_edge(u, v);
        return (long long) a * n + b;
    };

    for (int step = 0; step < steps; ++step) {
        double p = op(rng);
        if (p < 0.40) {
            // weight update
            int v = vd(rng);
            uint32_t nw = (uint32_t) wd(rng);
            w[v] = nw;
            tree.UpdateWeight(v, nw);
        } else if (p < 0.70) {
            // link attempt across components
            int u = vd(rng), v = vd(rng);
            if (u == v) continue;
            if (ref_connected(adj, u, v)) continue;
            parlay::sequence<std::pair<int, int>> link = {{u, v}};
            tree.BatchLink(link);
            adj[u].insert(v);
            adj[v].insert(u);
            edge_set.insert(edge_key(u, v));
        } else if (!edge_set.empty()) {
            // cut existing edge
            int idx = std::uniform_int_distribution<int>(0, (int) edge_set.size() - 1)(rng);
            auto it = edge_set.begin();
            std::advance(it, idx);
            long long key = *it;
            int u = (int) (key / n);
            int v = (int) (key % n);
            parlay::sequence<std::pair<int, int>> cut = {{u, v}};
            tree.BatchCut(cut);
            adj[u].erase(v);
            adj[v].erase(u);
            edge_set.erase(it);
        }

        // query check
        int a = vd(rng), b = vd(rng);
        bool conn = ref_connected(adj, a, b);
        EXPECT_EQ(tree.connected(a, b), conn);
        if (conn) {
            uint32_t expected = ref_path_sum(adj, w, a, b);
            uint32_t got = tree.PathQuery(a, b);
            if (got != expected) {
                auto path = ref_path_vertices(adj, a, b);
                std::cerr << "Mismatch details: step=" << step << " query=(" << a << "," << b << ")\n";
                std::cerr << "Path: ";
                for (int x : path) std::cerr << x << " ";
                std::cerr << "\nWeights: ";
                for (int x : path) std::cerr << w[x] << " ";
                std::cerr << "\nLeafDegrees: ";
                for (int x : path) std::cerr << tree.leaves[x].get_degree() << " ";
                std::cerr << "\nDeg(a)=" << tree.leaves[a].get_degree() << " Deg(b)=" << tree.leaves[b].get_degree() << "\n";
            }
            EXPECT_EQ(got, expected) << "step=" << step << " query=(" << a << "," << b << ")";
        }
    }
}

TEST(IntSumParallelUFOTreePathQuerySuite, random_dynamic_forest_multi_seed_multi_size) {
    auto run_trial = [&](int n, int steps, int k, uint32_t seed) {
        ParallelUFOTree<> tree(n, k);
        std::vector<std::unordered_set<int>> adj(n);
        std::unordered_set<long long> edge_set;
        std::vector<uint32_t> w(n, 1);

        std::mt19937 rng(seed);
        std::uniform_int_distribution<int> vd(0, n - 1);
        std::uniform_int_distribution<int> wd(1, 2000);
        std::uniform_real_distribution<double> op(0.0, 1.0);

        auto edge_key = [&](int u, int v) {
            auto [a, b] = norm_edge(u, v);
            return (long long) a * n + b;
        };

        for (int step = 0; step < steps; ++step) {
            double p = op(rng);
            if (p < 0.40) {
                int v = vd(rng);
                uint32_t nw = (uint32_t) wd(rng);
                w[v] = nw;
                tree.UpdateWeight(v, nw);
            } else if (p < 0.70) {
                int u = vd(rng), v = vd(rng);
                if (u == v) continue;
                if (ref_connected(adj, u, v)) continue;
                parlay::sequence<std::pair<int, int>> link = {{u, v}};
                tree.BatchLink(link);
                adj[u].insert(v);
                adj[v].insert(u);
                edge_set.insert(edge_key(u, v));
            } else if (!edge_set.empty()) {
                int idx = std::uniform_int_distribution<int>(0, (int) edge_set.size() - 1)(rng);
                auto it = edge_set.begin();
                std::advance(it, idx);
                long long key = *it;
                int u = (int) (key / n);
                int v = (int) (key % n);
                parlay::sequence<std::pair<int, int>> cut = {{u, v}};
                tree.BatchCut(cut);
                adj[u].erase(v);
                adj[v].erase(u);
                edge_set.erase(it);
            }

            int a = vd(rng), b = vd(rng);
            bool conn = ref_connected(adj, a, b);
            ASSERT_EQ(tree.connected(a, b), conn) << "n=" << n << " seed=" << seed << " step=" << step;
            if (conn) {
                uint32_t expected = ref_path_sum(adj, w, a, b);
                uint32_t got = tree.PathQuery(a, b);
                if (got != expected) {
                    auto path = ref_path_vertices(adj, a, b);
                    std::cerr << "Stress mismatch details: n=" << n << " seed=" << seed
                              << " step=" << step << " query=(" << a << "," << b << ")\n";
                    std::cerr << "Path: ";
                    for (int x : path) std::cerr << x << " ";
                    std::cerr << "\nWeights: ";
                    for (int x : path) std::cerr << w[x] << " ";
                    std::cerr << "\nLeafDegrees: ";
                    for (int x : path) std::cerr << tree.leaves[x].get_degree() << " ";
                    std::cerr << "\nDeg(a)=" << tree.leaves[a].get_degree()
                              << " Deg(b)=" << tree.leaves[b].get_degree() << "\n";
                }
                ASSERT_EQ(got, expected)
                    << "n=" << n << " seed=" << seed << " step=" << step
                    << " query=(" << a << "," << b << ")";
            }
        }
    };

    const std::vector<int> sizes = {12, 24, 40, 64};
    const std::vector<uint32_t> seeds = {1u, 42u, 777u, 20240528u, 987654321u};
    for (int n : sizes) {
        int k = std::max(4, n / 3);
        int steps = 12 * n;
        for (uint32_t seed : seeds) run_trial(n, steps, k, seed);
    }
}


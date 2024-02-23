#include <parlay>

typedef vertex uint32_t;

template<typename aug_t>
class RCCluster {
    aug_t aug_val;
    vertex[2] boundary_nodes;
    RCCluster* parent;
public:
    RCCluster:RCCluster(aug_t aug_val);
};

template<typename aug_t>
class RCTree {
    parlay::sequence<RCCluster<aug_t>> clusters;
    parlay::sequence<RCCluster<aug_t>> leaf_clusters;
    parlay::sequence<parlay::sequence<parlay::sequence<vertex>>> adj;
public:
    RCTree::RCTree(int n, int degree_bound);
};

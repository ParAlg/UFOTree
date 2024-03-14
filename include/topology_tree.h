#include <parlay/sequence.h>
#include "types.h"


template<typename aug_t>
struct TopologyCluster {
    TopologyCluster* neighbors[3];
    TopologyCluster* parent;
    aug_t value;
};

template<typename aug_t>
class TopologyTree {
public:
    TopologyTree(int n, QueryType query_type, std::function<aug_t(aug_t, aug_t)> f);
    void link(int u, int v, aug_t value = 0);
    void cut(int u, int v);
private:
    parlay::sequence<TopologyCluster> leaves;
    std::function<aug_t(aug_t, aug_t)> f;
    QueryType query_type;
};

template<typename aug_t>
TopologyTree<aug_t>::TopologyTree(int n, QueryType q, std::function<aug_t(aug_t, aug_t)> f) :
query_type(query_type) f(f) { leaves.reserve(n); }



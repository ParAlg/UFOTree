#include <vector>
#include "types.h"


template<typename DynamicTree, typename aug_t>
class TernarizedTree {
public:
    // Ternarized tree interface
    TernarizedTree(vertex_t n, QueryType q = PATH, std::function<aug_t(aug_t, aug_t)> f = [](aug_t x, aug_t y)->aug_t{return x + y;}, aug_t id = 0, aug_t dval = 0);
    void link(vertex_t u, vertex_t v, aug_t value = 1);
    void cut(vertex_t u, vertex_t v);
    bool connected(vertex_t u, vertex_t v);
private:
    // Underlying dynamic tree data structure
    DynamicTree tree;
    // Ternarization book-keeping
    vertex_t n;
    std::vector<vertex_t> free_vertex_ids;
    std::vector<vertex_t> vertex_id_map;
};

template<typename DynamicTree, typename aug_t>
TernarizedTree<DynamicTree, aug_t>::TernarizedTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d) : 
    tree(2*n, q, f, id, d), n(n) {}

template<typename DynamicTree, typename aug_t>
void TernarizedTree<DynamicTree, aug_t>::link(vertex_t u, vertex_t v, aug_t value) {
    // TODO: replace this with appropriate logic
    tree.link(u, v, value);
}

template<typename DynamicTree, typename aug_t>
void TernarizedTree<DynamicTree, aug_t>::cut(vertex_t u, vertex_t v) {
    // TODO: replace this with appropriate logic
    tree.cut(u, v);
}

template<typename DynamicTree, typename aug_t>
bool TernarizedTree<DynamicTree, aug_t>::connected(vertex_t u, vertex_t v) {
    return tree.connected(u,v);
}
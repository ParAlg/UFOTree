#include "../baselines/dynamic_trees/toptree-c-example/tree.c"
#include "../baselines/dynamic_trees/toptree-c-example/top_tree.c"
#include "../util/util.h"
#include "../util/types.h"

template <typename aug_t>
class TopTree {
    public:
        struct tree t;
        TopTree(size_t n, QueryType q, std::function<aug_t(aug_t,aug_t)> f, aug_t id, aug_t d_val);
        void link(vertex_t u, vertex_t v, aug_t weight);
        void cut(vertex_t u, vertex_t v);
};

template <typename aug_t>
TopTree<aug_t>::TopTree(size_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d_val){
    t = create_tree(n);
}

template <typename aug_t>
void TopTree<aug_t>::link(vertex_t u, vertex_t v, aug_t weight){
   link(&t.vertices[u], &t.vertices[v], weight);
}

//void TopTree<aug_t>::cut

// Given 2 vertices, call expose - traverse down the tree for each tt_node that has u and v as boundary vertices. Keep this going until tt_node->leaf = true.
template <typename aug_t>
void TopTree<aug_t>::cut(vertex_t u, vertex_t v){
    
}
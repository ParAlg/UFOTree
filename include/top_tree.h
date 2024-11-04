#include "../baselines/dynamic_trees/toptree-c-example/tree.c"
#include "../baselines/dynamic_trees/toptree-c-example/top_tree.c"
#include "../util/util.h"
#include "../util/types.h"

template <typename aug_t>
class TopTree {
    public:
        struct tree t;
        TopTree(size_t n, QueryType q = PATH, std::function<aug_t(aug_t,aug_t)> f = [] (int a, int b){return std::min(a,b);}, 
                aug_t id = std::numeric_limits<int>::max(), aug_t d_val = std::numeric_limits<int>::max());
        void link(vertex_t u, vertex_t v, aug_t weight);
        bool connected(vertex_t u, vertex_t v);
        void cut(vertex_t u, vertex_t v);
};

template <typename aug_t>
TopTree<aug_t>::TopTree(size_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d_val){
    t = create_tree(n);
}

template <typename aug_t>
bool TopTree<aug_t>::connected(vertex_t u, vertex_t v){
    tt_node* root1 = expose(&t.vertices[u]);
    tt_node* root2 = expose(&t.vertices[v]);
    auto to_ret = root1 && root1 == root2;
    deexpose(&t.vertices[u]);
    deexpose(&t.vertices[v]);
    return to_ret;
}
template <typename aug_t>
void TopTree<aug_t>::link(vertex_t u, vertex_t v, aug_t weight)
{
    struct vertex *ui = &t.vertices[u];
    struct vertex *vi = &t.vertices[v];
    tt_link(ui, vi, weight);
}

// void TopTree<aug_t>::cut

// Given 2 vertices, call expose - traverse down the tree for each tt_node that has u and v as boundary vertices. Keep this going until tt_node->leaf = true.
template <typename aug_t>
void TopTree<aug_t>::cut(vertex_t u, vertex_t v)
{
}
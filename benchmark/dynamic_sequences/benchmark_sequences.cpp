#include <parlay/parallel.h>
#include "pam.h"
#include "parett/sequence/treap/treap.hpp"
#include "parett/sequence/skip_list/skip_list.hpp"
#include "parett/sequence/parallel_skip_list/augmented_skip_list.hpp"
#include "parett/sequence/splay_tree/splay_tree.hpp"


using key_type = uint32_t;

struct entry {
    using key_t = key_type;
    static inline bool comp(key_t a, key_t b) { return a < b;}
    static size_t hash(key_t k) { return parlay::hash32(k);}
};

using treap_sequence  = pam_set<entry,pam_treap<entry>>;
using treap_seq_ops = treap_sequence::Seq_Tree;
using treap_map_ops = treap_sequence::Tree;


int main(int argc, char** argv) {
    int n = 10000000;
    parlay::sequence<key_type> v = parlay::tabulate(n, [&] (key_type i) { return i; });
    parlay::sequence<std::pair<int,int>> links = parlay::tabulate(n-1, [&] (size_t i) -> std::pair<int,int> { return {i,i+1}; });
    links = parlay::random_shuffle(links);
    pbbs::random randomness;
    parlay::internal::timer t;

    treap_sequence s1(v);
    treap_sequence::node* root = s1.root;
    t.start();
    for (int i = 0; i < n; i++) {
        auto split_result = treap_map_ops::split(root, v[i]);
        root = split_result.second;
    }
    t.next("PAM Treap Total Split Time");
    std::cout << std::endl;

    std::vector<treap::Node> treap;
    std::vector<treap::Node*> roots(n);
    for (int i = 0; i < n; i++)
        treap.emplace_back();
    t.start();
    for (int i = 0; i < links.size(); i++)
        treap::Node::Join(&treap[links[i].first], &treap[links[i].second]);
    t.next("ParETT Treap Total Join Time");
    t.start();
    for (int i = 0; i < n; i++)
        roots[i] = treap[i].GetRoot();
    t.next("ParETT Treap Total GetRoot Time");
    t.start();
    for (int i = 0; i < links.size(); i++)
        treap[links[i].first].Split();
    t.next("ParETT Treap Total Split Time");
    std::cout << std::endl;

    skip_list::Element* verts = pbbs::new_array_no_init<skip_list::Element>(n);
    std::vector <skip_list::Element*> sroots(n);
    for (int i = 0; i < n; i++)
        new (&verts[i]) skip_list::Element(randomness.ith_rand(i));
    t.start();
    for (int i = 0; i < links.size(); i++)
        skip_list::Element::Join(&verts[links[i].first], &verts[links[i].second]);
    t.next("ParETT SkipList Total Join Time");
    t.start();
    for (int i = 0; i < n; i++)
        sroots[i] = verts[i].FindRepresentative();
    t.next("ParETT SkipList Total GetRoot Time");
    t.start();
    for (int i = 0; i < links.size(); i++)
        verts[links[i].first].Split();
    t.next("ParETT SkipList Total Split Time");
    std::cout << std::endl;

    using Element = parallel_skip_list::AugmentedElement<int>;
    Element::aggregate_function = [&] (int x, int y) { return x+y; };
    Element* pverts = pbbs::new_array_no_init<Element>(n);
    std::vector <Element*> proots(n);
    for (int i = 0; i < n; i++)
        new (&pverts[i]) Element(randomness.ith_rand(i));
    t.start();
    for (int i = 0; i < links.size(); i++)
        Element::SequentialJoin2(&pverts[links[i].first], &pverts[links[i].second]);
    Element::SequentialJoin2(&pverts[n-1], &pverts[0]);
    t.next("ParETT ParSkipList Total Join Time");
    t.start();
    for (int i = 0; i < n; i++)
        proots[i] = pverts[i].FindRepresentative2();
    t.next("ParETT ParSkipList Total GetRoot Time");
    t.start();
    for (int i = 0; i < links.size(); i++)
        pverts[links[i].first].SequentialSplitRight();
    t.next("ParETT ParSkipList Total Split Time");
    std::cout << std::endl;

    std::vector<splay_tree::Node> splay_tree;
    for (int i = 0; i < n; i++)
        splay_tree.emplace_back();
    t.start();
    for (int i = 0; i < links.size(); i++)
        splay_tree::Node::Join(&splay_tree[links[i].first], &splay_tree[links[i].second]);
    t.next("ParETT SplayTree Total Join Time");
    t.start();
    for (int i = 0; i < n-1; i++)
        splay_tree[links[i].first].Split();
    t.next("ParETT SplayTree Total Split Time");
    std::cout << std::endl;
}
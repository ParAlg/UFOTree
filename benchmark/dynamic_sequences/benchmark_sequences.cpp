#include <parlay/parallel.h>
#include "pam.h"
#include "parett/sequence/treap/treap.hpp"


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
    parlay::internal::timer t;

    treap_sequence s1(v);
    treap_sequence::node* root = s1.root;
    t.start();
    for (int i = 0; i < n; i++) {
        auto split_result = treap_map_ops::split(root, v[i]);
        root = split_result.second;
    }
    t.next("PAM Treap Total Split Time");

    std::vector<treap::Node> treap;
    for (int i = 0; i < n; i++)
        treap.emplace_back();
    t.start();
    for (int i = 0; i < n-1; i++)
        treap::Node::Join(&treap[i], &treap[i+1]);
    t.next("ParETT Treap Total Join Time");
    t.start();
    for (int i = 0; i < n-1; i++)
        treap[i].Split();
    t.next("ParETT Treap Total Split Time");
}
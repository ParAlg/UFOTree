#include <parlay/parallel.h>
#include "pam.h"
// #include "baselines/sequence/parallel_treap/include/treap.hpp"


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
    size_t n = 1000000;
    parlay::sequence<key_type> v = parlay::tabulate(n, [&] (key_type i) { return i; });

    treap_sequence s1(v);
    treap_sequence::node* root = s1.root;
    parlay::internal::timer t;
    t.start();
    for (size_t i = 0; i < n; i++) {
        auto split_result = treap_map_ops::split(root, v[i]);
        root = split_result.second;
    }
    t.next("PAM Treap Total Split Time");


}
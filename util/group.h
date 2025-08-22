#pragma once
#include <parlay/parallel.h>
#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include "semisort.h"
#include "samplesort.h"
#include "hwy/base.h"


namespace dgbs {

template<typename C>
inline hwy::uint128_t make_hwy(C x, C y) {
    return (hwy::uint128_t) {(uintptr_t) y, (uintptr_t) x};
}

template<typename C>
inline std::pair<C, C> unpack_hwy(hwy::uint128_t x) {
    return std::make_pair((C) x.hi, (C) x.lo);
}

template <typename T>
auto group_by_key_inplace(parlay::sequence<T>& seq) {
    sample_sort_inplace(seq);
    auto starts = parlay::delayed_tabulate(seq.size(), [&] (int i) {
        if (i == 0 || seq[i].hi != seq[i-1].hi) return true;
        return false;
    });
    auto offsets = parlay::pack_index(starts);
    return parlay::tabulate(offsets.size(), [&] (size_t i) {
        size_t start = offsets[i];
        size_t end = i == offsets.size()-1 ? seq.size() : offsets[i+1];
        return std::make_pair(seq[offsets[i]].hi, parlay::make_slice(seq.begin() + start, seq.begin() + end));
    });
}

// template <typename K, typename V>
// auto group_by_key_inplace_old(parlay::sequence<std::pair<K, V>>& seq) {
//     seq = parlay::sort(seq);
//     auto starts = parlay::delayed_tabulate(seq.size(), [&] (int i) {
//         if (i == 0 || seq[i].first != seq[i-1].first) return true;
//         return false;
//     });
//     auto offsets = parlay::pack_index(starts);
//     return parlay::tabulate(offsets.size(), [&] (size_t i) {
//         size_t start = offsets[i];
//         size_t end = i == offsets.size()-1 ? seq.size() : offsets[i+1];
//         return std::make_pair(seq[offsets[i]].first, parlay::make_slice(seq.begin() + start, seq.begin() + end));
//     });
// }

template <typename K, typename V>
auto integer_group_by_key_inplace(parlay::sequence<std::pair<K, V>>& seq) {
    // parlay::integer_sort_inplace(seq, [&] (auto pair) { return (uintptr_t) pair.first; });
    auto seq_slice = parlay::make_slice(seq.begin(), seq.end());
    parlay::semisort_equal_inplace(seq_slice,
        [&] (auto x) { return x.first; },
        [&] (auto x) { return (uintptr_t) x; },
        [&] (auto x, auto y) { return x == y; }
    );
    auto starts = parlay::delayed_tabulate(seq.size(), [&] (int i) {
        if (i == 0 || seq[i].first != seq[i-1].first) return true;
        return false;
    });
    auto offsets = parlay::pack_index(starts);
    return parlay::tabulate(offsets.size(), [&] (size_t i) {
        size_t start = offsets[i];
        size_t end = i == offsets.size()-1 ? seq.size() : offsets[i+1];
        return std::make_pair(seq[offsets[i]].first, parlay::make_slice(seq.begin() + start, seq.begin() + end));
    });
}

}

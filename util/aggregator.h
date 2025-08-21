#pragma once
#include <utility>
#include <cstdint>
#include <parlay/parallel.h>
#include <parlay/sequence.h>


namespace dgbs {

template <class _T>
struct Aggregator {
    using T = _T;
    using TList = std::vector<T>;
    size_t padding;
    size_t num_workers;
    std::vector<TList> lists;
    Aggregator(size_t reserve_space = 0) {
        padding = (128 + sizeof(TList) - 1) / sizeof(TList);
        num_workers = parlay::num_workers();
        lists = std::vector<TList>(padding*num_workers);
        if (reserve_space > 0) {
            for (size_t i = 0; i < num_workers; ++i) {
                lists[padding*i].reserve(reserve_space);
            }
        }
    }

    void push_back(T x) {
        lists[padding*parlay::worker_id()].push_back(x);
    }

    bool empty() {
        for (size_t i = 0; i < num_workers; ++i)
            if (lists[padding*i].size() > 0)
                return false;
        return true;
    }

    template <class F>
    void for_all(F f) {
        parlay::parallel_for(0, num_workers, [&] (size_t i) {
            // for (size_t j = 0; j < lists[padding*i].size(); ++j) {
            //     f(lists[padding*i][j]);
            // }
            // parlay::parallel_for(0, lists[padding*i].size(), [&] (size_t j) {
            //   f(lists[padding*i][j]);
            // });
            // parlay::blocked_for(0, lists[padding*i].size(), 64, [&] (size_t _, size_t s, size_t e) {
            //   for (size_t j = s; j < e; j++) {
            //     f(lists[padding*i][j]);
            //   }
            // });
            if (lists[padding*i].size() < 64) {
              for (size_t j = 0; j < lists[padding*i].size(); j++) {
                f(lists[padding*i][j]);
              }
            } else {
              parlay::parallel_for(0, lists[padding*i].size(), [&] (size_t j) {
                f(lists[padding*i][j]);
              });
            }
        }, 1);
    }

    void clear() {
        parlay::parallel_for(0, num_workers, [&] (size_t i) {
            lists[padding*i].clear();
        }, 1);
    }

    parlay::sequence<T> to_sequence() {
        auto non_padding_lists = parlay::delayed_tabulate(num_workers, [&] (size_t i) {
            return lists[padding*i];
        });
        return parlay::flatten(non_padding_lists);
    }
};

}

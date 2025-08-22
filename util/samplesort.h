#include <algorithm>
#include <functional>
#include <random>

#include <parlay/parallel.h>
#include <parlay/primitives.h>
#include <parlay/random.h>
#include <parlay/slice.h>
#include <parlay/utilities.h>

#include <hwy/contrib/sort/vqsort-inl.h>
#include <hwy/contrib/sort/vqsort.h>
#include <hwy/highway.h>

// An efficient search tree for replacing binary-search on a sorted sequence.
// Returns the rank of the first element greater than or equal to the key.
// Reorganizes the values into a "heap ordering":
// i.e. with root at 0 and the children of position i are at 2*i+1 and 2*i+2.
// Significantly more efficient than binary search when tree fits in cache
// since it avoids conditionals.
// The number of pivots must be 2^i-1 (i.e. fully balanced tree)
template <typename T>
struct heap_tree {
 private:
  const long size;
  parlay::sequence<T> tree;
  const long levels;

  // converts from sorted sequence into a "heap indexed" tree
  void to_tree(const parlay::sequence<T>& In,
               long root, long l, long r) {
    size_t n = r - l;
    size_t m = l + n / 2;
    tree[root] = In[m];
    if (n == 1) return;
    to_tree(In, 2 * root + 1, l, m);
    to_tree(In, 2 * root + 2, m + 1, r);
  }
 public:
  // constructor
  heap_tree(const parlay::sequence<T>& keys) :
      size(keys.size()), tree(parlay::sequence<T>(size)),
      levels(parlay::log2_up(keys.size()+1)-1) {
    to_tree(keys, 0, 0, size);}

  // finds a key in the "heap indexed" tree
  template <typename Compare>
  inline int find(const T& key, Compare&& compare) {
    long j = 0;
    for (int k = 0; k < levels; k++) {
      j = 1 + 2 * j + compare(tree[j],key);
    }
    bool lt = compare(key,tree[j]);
    bool gt = compare(tree[j], key);
    j = 1 + 2 * j + !lt + -1*(!lt && !gt);
    return j - size;
  }

  // finds a key in the "heap indexed" tree
  inline int find(const T& key) {
    return find(key, std::less<>{});
  }
};


// **************************************************************
// Sample sort
// A generalization of quicksort to many pivots.
// This code picks up to 256 pivots by randomly selecting and
// then sorting them.
// It then puts the keys into buckets depending on which pivots
// they fall between and recursively sorts within the buckets.
// Makes use of a parlaylib bucket sort for the bucketing,
// and std::sort for the base case and sorting the pivots.
// **************************************************************

template <class Iter>
void VQSort(Iter b, Iter e) {
  hwy::VQSort(b, e-b, hwy::SortAscending());
}

template <typename Range, typename Less>
void sample_sort_(Range in, Range out, Less less, int level=1) {
  long n = in.size();

  // for the base case (small or recursion level greater than 2) use std::sort
  long cutoff = 2048;
  if (n <= cutoff || level > 2) {
    parlay::copy(in, out);
    // std::sort(out.begin(), out.end());
    VQSort(out.begin(), out.end());
    return;
  }

  // number of bits in bucket count (e.g. 8 would mean 256 buckets)
  int bits = std::min<long>(8, parlay::log2_up(n)-parlay::log2_up(cutoff)+1);
  long num_buckets = 1 << bits;

  // over-sampling ratio: keeps the buckets more balanced
  int over_ratio = 8;

  // create an over sample and sort it using std:sort
  parlay::random_generator gen;
  std::uniform_int_distribution<long> dis(0, n-1);
  auto oversample = parlay::tabulate(num_buckets * over_ratio, [&] (long i) {
    auto r = gen[i];
    return in[dis(r)];});
  std::sort(oversample.begin(), oversample.end());

  // sub sample to pick final pivots (num_buckets - 1 of them)
  auto pivots = parlay::tabulate(num_buckets-1, [&] (long i) {
    return oversample[(i+1)*over_ratio];});

  // put pivots into efficient search tree and find buckets id for the input keys
  heap_tree ss(pivots);
  auto bucket_ids = parlay::tabulate(n, [&] (long i) -> unsigned char {
    return ss.find(in[i], less);});

  // sort into the buckets
  auto [keys,offsets] = parlay::internal::count_sort(in, bucket_ids, num_buckets);

  // now recursively sort each bucket
  parlay::parallel_for(0, num_buckets, [&, &keys = keys, &offsets = offsets] (long i) {
    long first = offsets[i];
    long last = offsets[i+1];
    if (last-first == 0) return; // empty
    // are all equal, then can copy and quit
    if (i == 0 || i == num_buckets - 1 || less(pivots[i - 1], pivots[i]))
      sample_sort_(keys.cut(first,last), out.cut(first,last), less, level+1);
    else parlay::copy(keys.cut(first,last), out.cut(first,last));
  }, 1);
}

//  A wraper that calls sample_sort_
template <typename Range, typename Less = std::less<>>
auto sample_sort(Range& in, Less less = {}) {
  auto in_slice = parlay::make_slice(in.begin(), in.end());
  auto out = parlay::sequence<hwy::uint128_t>::uninitialized(in.size());
  auto out_slice = parlay::make_slice(out.begin(), out.end());
  sample_sort_(in_slice, out_slice, less);
  return out;
}

template <typename Range, typename Less = std::less<>>
auto sample_sort_inplace(Range& in, Less less = {}) {
  auto in_slice = parlay::make_slice(in.begin(), in.end());
  sample_sort_(in_slice, in_slice, less);
  return in;
}

#pragma once

#include <algorithm>
#include <unordered_map>
#include <vector>
#include <absl/container/flat_hash_map.h>
#include <parett/sequence/skip_list/skip_list.hpp>
#include <parett/utilities/hash_pair.hpp>
#include <parett/utilities/random.h>

using Element = skip_list::Element;
using std::pair;


namespace skip_list_ett {

class EulerTourTree {
 public:
  EulerTourTree(int _num_verts);
  // Note: this calls [finish()] on a [list_allocator<ETElement>]. Since
  // [list_allocator<>] functions are all static, this may cause errors if any
  // other [list_allocator<ETElement>]s are around in use (e.g. by other
  // [EulerTourTree] structures)
  // TODO(ttseng) deal with aforementioned issue better
  ~EulerTourTree();
  // Constructors and operators for C++11's rule of 5
  // TODO(ttseng) unimplemented
  EulerTourTree(const EulerTourTree&) = delete;
  EulerTourTree(EulerTourTree&&) = delete;
  EulerTourTree& operator=(const EulerTourTree&) = delete;
  EulerTourTree& operator=(EulerTourTree&&) = delete;

  bool IsConnected(int u, int v);
  void Link(int u, int v);
  void Cut(int u, int v);

  bool connected(int u, int v){return IsConnected(u,v);}
  void link(int u, int v){Link(u,v);}
  void cut(int u, int v){ Cut(u,v);}

  bool* BatchConnected(std::pair<int, int>* queries, int len);
  // Inserting all links in [links] must keep the graph acylic.
  void BatchLink(std::pair<int, int>* links, int len);
  // All edges in [cuts] must be in the graph, and no edges may be repeated.
  void BatchCut(std::pair<int, int>* cuts, int len);

  // Return the space being used by the skip list.
  size_t space();
  size_t count_nodes();
  size_t get_height();

 private:
  int num_verts;
  skip_list::Element* verts;
  absl::flat_hash_map<std::pair<int, int>, skip_list::Element*,
    HashIntPairStruct> edges;
  std::vector<skip_list::Element*> node_pool;
};


EulerTourTree::EulerTourTree(int _num_verts) : num_verts(_num_verts) {
  pbbs::random randomness;
  verts = pbbs::new_array_no_init<Element>(num_verts);
  for (int i = 0; i < num_verts; i++) {
    new (&verts[i]) Element(randomness.ith_rand(i));
    Element::Join(&verts[i], &verts[i]);
  }
  randomness = randomness.next();
  const int max_edges = 2 * (num_verts - 1);
  edges.reserve(max_edges);
  node_pool.reserve(max_edges);
  for (int i = 0; i < max_edges; i++) {
    node_pool.push_back(new Element(randomness.ith_rand(i)));
  }
}

EulerTourTree::~EulerTourTree() {
  for (auto it : edges) {
    delete it.second;
  }
  for (auto it : node_pool) {
    delete it;
  }
  pbbs::delete_array(verts, num_verts);
}

size_t EulerTourTree::space() {
  size_t max_space = sizeof(EulerTourTree);
  max_space += num_verts * sizeof(skip_list::Element);
  max_space += edges.size() * (sizeof(std::pair<int, int>) + sizeof(Element*)) // Size of key value pairs
    + edges.bucket_count() * (sizeof(void*) + sizeof(size_t)); // Space used by linked list
  for (auto element : edges) {
    max_space += element.second->calculate_size();
  }
  for (auto curr : node_pool) {
    max_space += (sizeof(skip_list::Element*) + curr->calculate_size());
  }
  return max_space;
}

size_t EulerTourTree::count_nodes() {
  size_t node_count = 0;
  for (auto element : edges)
    node_count += element.second->height;
  return node_count;
}

size_t EulerTourTree::get_height() {
  size_t max_height = 0;
  for (auto element : edges)
    if (element.second->height > max_height)
      max_height = element.second->height;
  return max_height;
}

bool EulerTourTree::IsConnected(int u, int v) {
  return verts[u].FindRepresentative() == verts[v].FindRepresentative();
}

void EulerTourTree::Link(int u, int v) {
  Element* uv = node_pool.back();
  node_pool.pop_back();
  Element* vu = node_pool.back();
  node_pool.pop_back();
  edges[std::make_pair(u, v)] = uv;
  edges[std::make_pair(v, u)] = vu;
  Element* u_left = &verts[u];
  Element* v_left = &verts[v];
  Element* u_right = u_left->Split();
  Element* v_right = v_left->Split();
  Element::Join(u_left, uv);
  Element::Join(uv, v_right);
  Element::Join(v_left, vu);
  Element::Join(vu, u_right);
}

void EulerTourTree::Cut(int u, int v) {
  auto uv_it = edges.find(std::make_pair(u, v));
  auto vu_it = edges.find(std::make_pair(v, u));
  Element* uv = uv_it->second;
  Element* vu = vu_it->second;
  edges.erase(uv_it);
  edges.erase(vu_it);
  Element* u_left = uv->GetPreviousElement();
  Element* v_left = vu->GetPreviousElement();
  Element* v_right = uv->Split();
  Element* u_right = vu->Split();
  u_left->Split();
  v_left->Split();
  Element::Join(u_left, u_right);
  Element::Join(v_left, v_right);
  node_pool.push_back(uv);
  node_pool.push_back(vu);
}

bool* EulerTourTree::BatchConnected(pair<int, int>* queries, int len) {
  bool* ans = pbbs::new_array_no_init<bool>(len);
  for (int i = 0; i < len; i++) {
    ans[i] = IsConnected(queries[i].first, queries[i].second);
  }
  return ans;
}

void EulerTourTree::BatchLink(pair<int, int>* links, int len) {
  for (int i = 0; i < len; i++) {
    Link(links[i].first, links[i].second);
  }
}

void EulerTourTree::BatchCut(pair<int, int>* cuts, int len) {
  for (int i = 0; i < len; i++) {
    Cut(cuts[i].first, cuts[i].second);
  }
}

} // namespace skip_list_ett

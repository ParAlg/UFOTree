#include <dynamic_trees/link_cut_tree/include/link_cut_tree.hpp>
#include <absl/container/flat_hash_map.h>

namespace link_cut_tree {

class Node {
 public:
  Node();
  Node(Node* _par, Node* left, Node *right);

  Node* get_root();
  void cut(Node* neighbor);
  void cut_from_par();
  void link(Node* child);
  Node* lca(Node* other);
  void evert(); // reroot
  int path_to_root_query();

 private:
  Node* par; // parent
  Node* c[2]; // children
  bool flip; // whether children are reversed; used for evert()

  Node* get_real_par();
  void rot();
  void splay();
  Node* expose();
  void fix_c();
  void push_flip();
};

Node::Node(Node* _par, Node* left, Node *right)
  : par(_par), c{left, right}, flip(0) {
  fix_c();
}

Node::Node() : Node(nullptr, nullptr, nullptr) {}

Node* Node::get_real_par() {
  return par != nullptr && this != par->c[0] && this != par->c[1] ? nullptr : par;
}

void Node::fix_c() {
  for (int i = 0; i < 2; i++)
    if (c[i] != nullptr)
      c[i]->par = this;
}

void Node::push_flip() {
  if (flip) {
    flip = 0;
    std::swap(c[0], c[1]);
    for (int i = 0; i < 2; i++)
      if (c[i] != nullptr)
        c[i]->flip ^= 1;
  }
}

void Node::rot() { // rotate v towards its parent; v must have real parent
  Node* p = get_real_par();
  par = p->par;
  if (par != nullptr)
    for (int i = 0; i < 2; i++)
      if (par->c[i] == p) {
        par->c[i] = this;
        par->fix_c();
      }
  const bool rot_dir = this == p->c[0];
  p->c[!rot_dir] = c[rot_dir];
  c[rot_dir] = p;
  p->fix_c();
  fix_c();
}

void Node::splay() {
  Node* p, * gp;
  push_flip(); // guarantee flip bit isn't set after calling splay()
  while ((p = get_real_par()) != nullptr) {
    gp = p->get_real_par();
    if (gp != nullptr)
      gp->push_flip();
    p->push_flip();
    push_flip();
    if (gp != nullptr)
      ((gp->c[0] == p) == (p->c[0] == this) ? p : this)->rot();
    rot();
  }
}

// returns 1st vertex encountered that was originally in same path as root (used
// for LCA)
Node* Node::expose() {
  Node* ret = this;
  for (Node* curr = this, * pref = nullptr; curr != nullptr;
       ret = curr, pref = this, curr = par) {
    curr->splay();
    curr->c[1] = pref;
    curr->fix_c();
    splay();
  }
  return ret;
}

void Node::evert() {
  expose();
  flip ^= 1;
  push_flip();
}

Node* Node::get_root() {
  expose();
  Node* root = this;
  push_flip();
  while (root->c[0] != nullptr) {
    root = root->c[0];
    root->push_flip();
  }
  root->splay();
  return root;
}

void Node::cut_from_par() {
  expose();
  c[0] = c[0]->par = nullptr;
  fix_c();
}

void Node::cut(Node* neighbor) {
  neighbor->evert();
  evert();
  neighbor->par = nullptr;
  for (int i = 0; i < 2; i++)
    if (c[i] == neighbor)
      c[i] = nullptr;
  fix_c();
}

void Node::link(Node* child) {
  child->evert();
  expose();
  child->par = this;
}

Node* Node::lca(Node* other) {
  expose();
  return other->expose();
}

LinkCutTree::LinkCutTree(int _num_verts) : num_verts(_num_verts) {
  verts = new Node[num_verts];
}

LinkCutTree::~LinkCutTree() {
  delete[] verts;
}

void LinkCutTree::link(vertex_t u, vertex_t v) {
  verts[u].link(&verts[v]);
}

void LinkCutTree::cut(vertex_t u, vertex_t v) {
  verts[u].cut(&verts[v]);
}

bool LinkCutTree::connected(vertex_t u, vertex_t v) {
  return verts[u].get_root() == verts[v].get_root();
}

bool* LinkCutTree::BatchConnected(std::pair<int, int>* queries, int len) {
  bool* ans = new bool[len];
  for (int i = 0; i < len; i++)
    ans[i] = verts[queries[i].first].get_root() == verts[queries[i].second].get_root();
  return ans;
}

void LinkCutTree::BatchLink(std::pair<int, int>* links, int len) {
  for (int i = 0; i < len; i++)
    verts[links[i].first].link(&verts[links[i].second]);
}

void LinkCutTree::BatchCut(std::pair<int, int>* cuts, int len) {
  for (int i = 0; i < len; i++)
    verts[cuts[i].first].cut(&verts[cuts[i].second]);
}


// ============= Link Cut Tree With Integer Path Queries Below This Point ===============

class NodeInt {
 public:
  NodeInt();

  void link(NodeInt* child, int weight);
  void cut(NodeInt* neighbor);
  void evert(); // reroot
  NodeInt* get_root();
  int path_query(NodeInt* other);

 private:
  NodeInt* par; // parent
  NodeInt* c[2]; // children
  NodeInt* n[2]; // store the nodes directly left and right in the splay tree by key
  int max; // maintain the maximum edge weight in the splay tree subtree rooted at this
  bool flip; // whether children are reversed; used for evert()
  absl::flat_hash_map<NodeInt*, int> edges; // store the weights of edges to all neighbors

  NodeInt* get_real_par();
  NodeInt* get_leftmost();
  void rot();
  void splay();
  NodeInt* expose();
  void fix_c();
  void recompute_max();
  void push_flip();
};

NodeInt::NodeInt() : par(nullptr), c{nullptr, nullptr}, n{nullptr, nullptr},
    max(MIN_INT), flip(0), edges() {}

NodeInt* NodeInt::get_real_par() {
  return par != nullptr && this != par->c[0] && this != par->c[1] ? nullptr : par;
}

NodeInt* NodeInt::get_leftmost() { // only called when this is root of splay tree
  NodeInt* left = this;
  push_flip();
  while (left->c[0] != nullptr) {
    left = left->c[0];
    left->push_flip();
  }
  return left;
}

void NodeInt::fix_c() {
  for (int i = 0; i < 2; i++)
    if (c[i] != nullptr)
      c[i]->par = this;
}

void NodeInt::recompute_max() {
  max = MIN_INT;
  if (n[0] != nullptr) max = std::max(max, edges[n[0]]);
  if (n[1] != nullptr) max = std::max(max, edges[n[1]]);
  for (int i = 0; i < 2; i++)
    if (c[i] != nullptr)
      max = std::max(max, c[i]->max);
}

void NodeInt::push_flip() {
  if (flip) {
    flip = 0;
    std::swap(c[0], c[1]);
    std::swap(n[0], n[1]);
    for (int i = 0; i < 2; i++)
      if (c[i] != nullptr)
        c[i]->flip ^= 1;
  }
}

void NodeInt::rot() { // rotate v towards its parent; v must have real parent
  NodeInt* p = get_real_par();
  par = p->par;
  if (par != nullptr)
    for (int i = 0; i < 2; i++)
      if (par->c[i] == p) {
        par->c[i] = this;
        par->fix_c();
      }
  const bool rot_dir = this == p->c[0];
  p->c[!rot_dir] = c[rot_dir];
  c[rot_dir] = p;
  p->fix_c();
  p->recompute_max();
  fix_c();
  recompute_max();
}

void NodeInt::splay() {
  NodeInt* p, * gp;
  push_flip(); // guarantee flip bit isn't set after calling splay()
  while ((p = get_real_par()) != nullptr) {
    gp = p->get_real_par();
    if (gp != nullptr)
      gp->push_flip();
    p->push_flip();
    push_flip();
    if (gp != nullptr)
      ((gp->c[0] == p) == (p->c[0] == this) ? p : this)->rot();
    rot();
  }
}

// returns the root of the tree
NodeInt* NodeInt::expose() {
  NodeInt* curr = this;
  NodeInt* head = nullptr;
  while (curr) {
    curr->splay();
    curr->c[1] = head;
    if (curr->n[1]) {
      curr->n[1]->push_flip();
      curr->n[1]->n[0] = nullptr;
      curr->n[1]->recompute_max();
    }
    curr->n[1] = head;
    if (head) {
      head->n[0] = curr;
      head->recompute_max();
    }
    curr->fix_c();
    curr->recompute_max();
    head = curr->get_leftmost();
    head->splay();
    curr = head->par;
  }
  return head;
}

void NodeInt::evert() {
  NodeInt* head = expose();
  head->flip ^= 1;
  head->push_flip();
}

NodeInt* NodeInt::get_root() {
  return expose();
}

int NodeInt::path_query(NodeInt* other) {
  evert();
  other->expose();
  return max;
}

void NodeInt::cut(NodeInt* neighbor) {
  neighbor->evert();
  evert();
  neighbor->push_flip();
  push_flip();
  neighbor->c[0] = nullptr;
  neighbor->n[0] = nullptr;
  par = nullptr;
  n[1] = nullptr;
  this->edges.erase(neighbor);
  neighbor->edges.erase(this);
}

void NodeInt::link(NodeInt* child, int weight) {
  this->edges.insert({child, weight});
  child->edges.insert({this, weight});
  child->evert();
  child->splay();
  expose();
  child->par = this;
}

LinkCutTreeInt::LinkCutTreeInt(int _num_verts) : num_verts(_num_verts) {
  verts = new NodeInt[num_verts];
}

LinkCutTreeInt::~LinkCutTreeInt() {
  delete[] verts;
}

void LinkCutTreeInt::link(vertex_t u, vertex_t v, int weight) {
  verts[u].link(&verts[v], weight);
}

void LinkCutTreeInt::cut(vertex_t u, vertex_t v) {
  verts[u].cut(&verts[v]);
}

bool LinkCutTreeInt::connected(vertex_t u, vertex_t v) {
  bool conn = verts[u].get_root() == verts[v].get_root();
  return conn;
}

int LinkCutTreeInt::path_query(vertex_t u, vertex_t v) {
  return verts[u].path_query(&verts[v]);
}

} // namespace link_cut_tree

#pragma once
#include <stdint.h>


namespace dgbs {

typedef uint32_t vertex_t;
static vertex_t NONE = -1;

struct empty_t {
};
static empty_t empty;

typedef uint64_t edge_t;

enum UpdateType {
  INSERT,
  DELETE
};

struct Edge {
public:
  vertex_t src;
  vertex_t dst;

  bool operator==(const Edge& other) const {
    return src == other.src && dst == other.dst;
  }
};

struct Update {
  UpdateType type;
  Edge edge;
};

enum QueryType {
  CONNECTIVITY,
  PATH,
  SUBTREE
};

struct Query {
  vertex_t u;
  vertex_t v;
};

};

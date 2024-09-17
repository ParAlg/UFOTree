#pragma once
#include <stdint.h>


typedef uint32_t vertex_t;
static vertex_t NONE = -1;
static vertex_t MARK = 0;

struct empty_t {
};
static empty_t empty;

typedef uint64_t edge_t;

enum UpdateType {
  INSERT,
  DELETE
};

struct Edge {
  vertex_t src;
  vertex_t dst;
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

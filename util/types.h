#pragma once
#include <stdint.h>
#include <utility>


typedef uint32_t vertex_t;
static vertex_t NONE = -1;

typedef std::pair<vertex_t, vertex_t> edge_t;

struct empty_t {
};
static empty_t empty;

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

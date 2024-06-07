#pragma once
#include <stdint.h>


typedef uint32_t vertex_t;

enum QueryType {
  PATH,
  SUBTREE
};

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

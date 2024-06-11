#pragma once
#include<cstdint>


typedef uint32_t vertex_t;
#define MAX_VERTEX_T (std::numeric_limits<uint32_t>::max())

enum QueryType {
  PATH,
  SUBTREE
};

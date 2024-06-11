#pragma once
#include <math.h>
#include "types.h"

#define MAX_VERTEX_T (std::numeric_limits<uint32_t>::max())

static int max_tree_height(vertex_t n) {
    return ceil(log2(n) / log2(1.2));
}

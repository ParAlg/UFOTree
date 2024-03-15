#pragma once
#include <math.h>
#include "types.h"


static int max_tree_height(vertex_t n) {
    return ceil(log2(n) / log2(1.2));
}

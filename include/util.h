#pragma once
#include <math.h>
#include "types.h"


#define MAX_VERTEX_T (std::numeric_limits<uint32_t>::max())

static int max_tree_height(vertex_t n) {
    return ceil(log2(n) / log2(1.2));
}

#define CAS(OBJ, EXP, DES) std::atomic_compare_exchange_weak(OBJ, EXP, DES)

// #define START_TIMER(X) auto X = std::chrono::high_resolution_clock::now()
// #define STOP_TIMER(X, T) T += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now()-X).count()
// #define PRINT_TIMER(S, T) std::cout << "    " << S << " (ms): " << T/1000000 << std::endl

#define START_TIMER(X) ;
#define STOP_TIMER(X, T) ;
#define PRINT_TIMER(S, T) ;

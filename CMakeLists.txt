cmake_minimum_required(VERSION 3.15)
project(PracticalDynamicTrees)
include (FetchContent)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# add_compile_options(-fsanitize=address)
# add_link_options(-fsanitize=address)
# add_compile_options(-fsanitize=undefined)
# add_link_options(-fsanitize=undefined)

FetchContent_Declare(parlaylib
  GIT_REPOSITORY  https://github.com/cmuparlay/parlaylib.git
  GIT_TAG         master
)
FetchContent_GetProperties(parlaylib)
if(NOT parlaylib_POPULATED)
  FetchContent_Populate(parlaylib)  
  add_subdirectory(${parlaylib_SOURCE_DIR} EXCLUDE_FROM_ALL)
endif()

FetchContent_Declare(CPAM
  GIT_REPOSITORY  https://github.com/ParAlg/CPAM.git
  GIT_TAG         cmake
)
FetchContent_GetProperties(CPAM)
if(NOT CPAM_POPULATED)
  FetchContent_Populate(CPAM)  
  message(${cpam_SOURCE_DIR})
  add_subdirectory(${cpam_SOURCE_DIR} EXCLUDE_FROM_ALL)
endif()

FetchContent_Declare(googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
FetchContent_MakeAvailable(googletest)

FetchContent_Declare(abseil-cpp
  GIT_REPOSITORY  git@github.com:abseil/abseil-cpp.git
  GIT_TAG         master
)
FetchContent_MakeAvailable(abseil-cpp)

add_executable(dynamic_tree_tests
  test/test_runner.cpp
  test/link_cut_tree_tests.cpp
  test/top_tree_tests.cpp
  test/rc_tree_tests.cpp
  test/topology_tree_tests.cpp
  test/ufo_tree_tests.cpp
  test/ternarization_tests.cpp
  baselines/parett/dynamic_trees/link_cut_tree/src/link_cut_tree.cpp
)
target_link_libraries(dynamic_tree_tests PRIVATE parlay GTest::gtest_main absl::flat_hash_map)
target_include_directories(dynamic_tree_tests PRIVATE include util baselines)

# add_executable(parallel_dynamic_tree_tests
#   test/parallel_test_runner.cpp
#   test/parallel_topology_tree_tests.cpp
#   test/parallel_ufo_tree_tests.cpp
# )
# target_link_libraries(parallel_dynamic_tree_tests PRIVATE parlay CPAM aspen GTest::gtest_main)
# target_include_directories(parallel_dynamic_tree_tests PRIVATE include include/forest_structures util)

add_executable(benchmark_update_speed
  benchmark/dynamic_trees/benchmark_update_speed.cpp
  baselines/parett/dynamic_trees/link_cut_tree/src/link_cut_tree.cpp
  baselines/parett/dynamic_trees/euler_tour_tree/src/skip_list_ett.cpp
  baselines/parett/sequence/skip_list/src/augmented_skip_list.cpp
  baselines/parett/dynamic_trees/euler_tour_tree/src/splay_tree_ett.cpp
  baselines/parett/sequence/splay_tree/src/splay_tree.cpp
)
target_link_libraries(benchmark_update_speed PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_update_speed PRIVATE include util baselines)

add_executable(benchmark_query_speed
  benchmark/dynamic_trees/benchmark_query_speed.cpp
  baselines/parett/dynamic_trees/link_cut_tree/src/link_cut_tree.cpp
  baselines/parett/dynamic_trees/euler_tour_tree/src/skip_list_ett.cpp
  baselines/parett/sequence/skip_list/src/augmented_skip_list.cpp
  baselines/parett/dynamic_trees/euler_tour_tree/src/splay_tree_ett.cpp
  baselines/parett/sequence/splay_tree/src/splay_tree.cpp
)
target_link_libraries(benchmark_query_speed PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_query_speed PRIVATE include util baselines)

add_executable(benchmark_peak_space
  benchmark/dynamic_trees/benchmark_peak_space.cpp
  baselines/parett/dynamic_trees/link_cut_tree/src/link_cut_tree.cpp
  baselines/parett/dynamic_trees/euler_tour_tree/src/skip_list_ett.cpp
  baselines/parett/sequence/skip_list/src/augmented_skip_list.cpp
  baselines/parett/dynamic_trees/euler_tour_tree/src/splay_tree_ett.cpp
  baselines/parett/sequence/splay_tree/src/splay_tree.cpp
)
target_link_libraries(benchmark_peak_space PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_peak_space PRIVATE include util baselines)

add_executable(benchmark_graph_statistics
  benchmark/graph_algorithms/benchmark_graph_statistics.cpp
  baselines/parett/dynamic_trees/link_cut_tree/src/link_cut_tree.cpp
  baselines/parett/dynamic_trees/euler_tour_tree/src/skip_list_ett.cpp
  baselines/parett/sequence/skip_list/src/augmented_skip_list.cpp
  baselines/parett/dynamic_trees/euler_tour_tree/src/splay_tree_ett.cpp
  baselines/parett/sequence/splay_tree/src/splay_tree.cpp
)
target_link_libraries(benchmark_graph_statistics PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_graph_statistics PRIVATE include util baselines)

add_executable(benchmark_graph_update_speed
  benchmark/graph_algorithms/benchmark_graph_update_speed.cpp
  baselines/parett/dynamic_trees/link_cut_tree/src/link_cut_tree.cpp
  baselines/parett/dynamic_trees/euler_tour_tree/src/skip_list_ett.cpp
  baselines/parett/sequence/skip_list/src/augmented_skip_list.cpp
  baselines/parett/dynamic_trees/euler_tour_tree/src/splay_tree_ett.cpp
  baselines/parett/sequence/splay_tree/src/splay_tree.cpp
)
target_link_libraries(benchmark_graph_update_speed PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_graph_update_speed PRIVATE include util baselines)

add_executable(benchmark_graph_algorithms
  benchmark/graph_algorithms/benchmark_graph_algorithms.cpp
  baselines/parett/dynamic_trees/link_cut_tree/src/link_cut_tree.cpp
  baselines/parett/dynamic_trees/euler_tour_tree/src/skip_list_ett.cpp
  baselines/parett/sequence/skip_list/src/augmented_skip_list.cpp
  baselines/parett/dynamic_trees/euler_tour_tree/src/splay_tree_ett.cpp
  baselines/parett/sequence/splay_tree/src/splay_tree.cpp
)
target_link_libraries(benchmark_graph_algorithms PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_graph_algorithms PRIVATE include util baselines)

# add_executable(parallel_benchmark
#   benchmark/parallel_benchmark_runner.cpp
#   test/parallel_topology_tree_tests.cpp
#   test/parallel_ufo_tree_tests.cpp
#   baselines/dynamic_trees/parallel_euler_tour_tree/src/euler_tour_tree.cpp
#   baselines/sequence/parallel_skip_list/src/augmented_skip_list.cpp
# )
# target_link_libraries(parallel_benchmark PRIVATE parlay CPAM aspen GTest::gtest_main)
# target_include_directories(parallel_benchmark PRIVATE include include/forest_structures util baselines)

# add_executable(incremental_MST
#   test/incremental_mst_test.cpp
# )
# target_link_libraries(incremental_MST PRIVATE parlay absl::flat_hash_map)
# target_include_directories(incremental_MST PRIVATE include util)

add_executable(benchmark_sequences
  benchmark/dynamic_sequences/benchmark_sequences.cpp
  # baselines/sequence/parallel_treap/src/treap.cpp
)
target_link_libraries(benchmark_sequences PRIVATE parlay)
target_include_directories(benchmark_sequences PRIVATE include util baselines baselines/pam)

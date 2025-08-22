cmake_minimum_required(VERSION 3.15)
cmake_policy(SET CMP0135 NEW)
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
# add_compile_options(-fsanitize=thread)
# add_link_options(-fsanitize=thread)

FetchContent_Declare(parlaylib
  GIT_REPOSITORY  https://github.com/cmuparlay/parlaylib.git
  GIT_TAG         master
)
FetchContent_GetProperties(parlaylib)
if(NOT parlaylib_POPULATED)
  FetchContent_Populate(parlaylib)  
  add_subdirectory(${parlaylib_SOURCE_DIR} EXCLUDE_FROM_ALL)
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

FetchContent_Declare(parallel-ett
    GIT_REPOSITORY  git@github.com:quintendeman/parallel-euler-tour-tree.git
    GIT_TAG         master
)
FetchContent_MakeAvailable(parallel-ett)

FetchContent_Declare(pam
    GIT_REPOSITORY  git@github.com:happyrogue2307/PAM.git
    GIT_TAG         master
)
FetchContent_MakeAvailable(pam)

set(GTEST_ROOT ${CMAKE_CURRENT_BINARY_DIR}/_deps/googletest-src)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)

FetchContent_Declare(highway
  GIT_REPOSITORY  https://github.com/google/highway/
  GIT_TAG         master
)
FetchContent_MakeAvailable(highway)

add_library(dgbs INTERFACE)
include(GNUInstallDirs)
set(DGBS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(dgbs INTERFACE
  $<BUILD_INTERFACE:${DGBS_INCLUDE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_include_directories(dgbs INTERFACE include util baselines)

add_executable(dynamic_tree_tests
  test/test_runner.cpp
  test/euler_tour_tree_tests.cpp
  test/link_cut_tree_tests.cpp
  test/top_tree_tests.cpp
  test/rc_tree_tests.cpp
  test/topology_tree_tests.cpp
  test/ufo_tree_tests.cpp
  test/ternarization_tests.cpp
)
target_link_libraries(dynamic_tree_tests PRIVATE parlay GTest::gtest_main absl::flat_hash_map)
target_include_directories(dynamic_tree_tests PRIVATE include util baselines)

add_executable(benchmark_update_speed
  benchmark/dynamic_trees/benchmark_update_speed.cpp
)
target_link_libraries(benchmark_update_speed PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_update_speed PRIVATE include util baselines)

add_executable(benchmark_query_speed
  benchmark/dynamic_trees/benchmark_query_speed.cpp
)
target_link_libraries(benchmark_query_speed PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_query_speed PRIVATE include util baselines)

add_executable(benchmark_peak_space
  benchmark/dynamic_trees/benchmark_peak_space.cpp
)
target_link_libraries(benchmark_peak_space PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_peak_space PRIVATE include util baselines)

add_executable(benchmark_graph_statistics
  benchmark/graph_algorithms/benchmark_graph_statistics.cpp
)
target_link_libraries(benchmark_graph_statistics PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_graph_statistics PRIVATE include util baselines)

add_executable(benchmark_graph_update_speed
  benchmark/graph_algorithms/benchmark_graph_update_speed.cpp
)
target_link_libraries(benchmark_graph_update_speed PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_graph_update_speed PRIVATE include util baselines)

add_executable(benchmark_graph_algorithms
  benchmark/graph_algorithms/benchmark_graph_algorithms.cpp
)
target_link_libraries(benchmark_graph_algorithms PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_graph_algorithms PRIVATE include util baselines)

add_executable(parallel_dynamic_tree_tests
  test/test_runner.cpp
  test/parallel_mis_tests.cpp
  test/parallel_topology_tree_tests.cpp
  test/parallel_ufo_tree_tests.cpp
  test/par_rc_tests.cpp
)
target_link_libraries(parallel_dynamic_tree_tests PRIVATE parlay GTest::gtest_main pam hwy hwy_contrib)
target_include_directories(parallel_dynamic_tree_tests PRIVATE include include/forest_structures util)

add_executable(parallel_benchmark
  benchmark/parallel/parallel_benchmark_runner.cpp
)
target_link_libraries(parallel_benchmark PRIVATE parlay GTest::gtest_main par-ett pam hwy hwy_contrib)
target_include_directories(parallel_benchmark PRIVATE benchmark/parallel include include/forest_structures util)
target_compile_options(parallel_benchmark PRIVATE -O3 -march=native)

# add_executable(incremental_MST
#   test/incremental_mst_test.cpp
# )
# target_link_libraries(incremental_MST PRIVATE parlay absl::flat_hash_map)
# target_include_directories(incremental_MST PRIVATE include util)

# add_executable(benchmark_sequences
#   benchmark/dynamic_sequences/benchmark_sequences.cpp
# )
# target_link_libraries(benchmark_sequences PRIVATE parlay par-ett)
# target_include_directories(benchmark_sequences PRIVATE include util baselines baselines/pam)

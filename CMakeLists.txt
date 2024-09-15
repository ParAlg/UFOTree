cmake_minimum_required(VERSION 3.15)
project(PracticalDynamicTrees)
include (FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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

add_executable(dynamic_tree_tests
  test/test_runner.cpp
  test/rc_tree_tests.cpp
  test/topology_tree_tests.cpp
  test/ufo_tree_tests.cpp
  test/ternarization_tests.cpp
)
target_link_libraries(dynamic_tree_tests PRIVATE parlay GTest::gtest_main absl::flat_hash_map)
target_include_directories(dynamic_tree_tests PRIVATE include)

add_executable(benchmark_update_speed
  benchmark/benchmark_update_speed.cpp
  baselines/dynamic_trees/euler_tour_tree/src/skip_list_ett.cpp
  baselines/sequence/skip_list/src/augmented_skip_list.cpp
)
target_link_libraries(benchmark_update_speed PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_update_speed PRIVATE include baselines)

add_executable(benchmark_query_speed
  benchmark/benchmark_query_speed.cpp
  baselines/dynamic_trees/euler_tour_tree/src/skip_list_ett.cpp
  baselines/sequence/skip_list/src/augmented_skip_list.cpp
)
target_link_libraries(benchmark_query_speed PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_query_speed PRIVATE include baselines)

add_executable(benchmark_peak_space
  benchmark/benchmark_peak_space.cpp
  baselines/dynamic_trees/euler_tour_tree/src/skip_list_ett.cpp
  baselines/sequence/skip_list/src/augmented_skip_list.cpp
)
target_link_libraries(benchmark_peak_space PRIVATE parlay absl::flat_hash_map)
target_include_directories(benchmark_peak_space PRIVATE include baselines)

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

add_executable(dynamic_tree_tests
  test/test_runner.cpp
  test/topology_tree_tests.cpp
  test/ufo_tree_tests.cpp
  test/test_rc_tree.cpp
)
target_link_libraries(dynamic_tree_tests PRIVATE parlay GTest::gtest_main)

add_executable(benchmark_dynamic_trees
  benchmark/benchmark_runner.cpp
  baselines/dynamic_trees/euler_tour_tree/src/skip_list_ett.cpp
  baselines/sequence/skip_list/src/augmented_skip_list.cpp
)
target_link_libraries(benchmark_dynamic_trees PRIVATE parlay)
target_include_directories(benchmark_dynamic_trees PRIVATE baselines)

cmake_minimum_required(VERSION 3.26)
project(smooth_backfitting_library)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math")

include_directories(include)
include_directories(extern)
#target_link_directories(ex3 PRIVATE include)
add_executable(ex2 examples/example_2.cpp)
add_executable(ex4 example_4.cpp)
add_executable(ex3 examples/example_3.cpp)
add_executable(SBF_test tests/smooth_backfitting_core_test.cpp
        tests/partially_linear_SBF_test.cpp)
add_executable(PL_SBF tests/partially_linear_SBF_test.cpp)


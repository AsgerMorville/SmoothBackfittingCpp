cmake_minimum_required(VERSION 3.26)
project(smooth_backfitting_library)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math")

include_directories(include)
include_directories(extern)
#target_link_directories(ex3 PRIVATE include)


add_executable(ex1 examples/example_1.cpp)
add_executable(ex2 examples/example_2.cpp)


add_executable(add_func_test tests/additive_function_test.cpp)
add_executable(random_generator_test tests/random_generator_test.cpp)
add_executable(SBF_test tests/smooth_backfitting_core_test.cpp)
add_executable(PL_SBF tests/partially_linear_SBF_test.cpp)


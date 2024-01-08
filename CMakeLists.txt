cmake_minimum_required(VERSION 3.26)
project(smooth_backfitting_library)

set(CMAKE_CXX_STANDARD 20)

include_directories(extern)
include_directories(include)

add_executable(example_1 examples/example_1.cpp
        include/additive_function.h
        tests/additive_function_test.cpp
        tests/smooth_backfitting_core_test.cpp)
add_executable(example_3 examples/example_3.cpp)
add_executable(additive_function_test tests/additive_function_test.cpp
        include/quantile.h)
add_executable(smooth_backfitting_core_test tests/smooth_backfitting_core_test.cpp)


cmake_minimum_required(VERSION 3.15)

option(ENABLE_TESTS "" OFF)
option(ENABLE_EXAMPLES "" OFF)
option(ENABLE_EXERCISES "" OFF)
add_subdirectory(RAJA)

# set(ENABLE_OPENMP ON)
set(ENABLE_BENCHMARKS OFF)
set(ENABLE_TESTS OFF)
set(ENABLE_GMOCK OFF)
option(ENABLE_EXAMPLES "" OFF)
add_subdirectory(Umpire)

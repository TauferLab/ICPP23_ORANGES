cmake_minimum_required(VERSION 3.15)

project(FIDO)

set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "")
find_package(OpenMP REQUIRED)

add_subdirectory(ESSENS-master)
add_subdirectory(GRAAL)

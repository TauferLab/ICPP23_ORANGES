set(CMAKE_CXX_COMPILER "dpcpp")

cmake_minimum_required(VERSION 3.15)

project(FIDO)

set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

#set(TRACK_RUNTIME OFF CACHE BOOL "")

include(GNUInstallDirs)

#add_subdirectory(ESSENS-master)
add_subdirectory(submodules/ESSENS)
add_subdirectory(fido)

install(EXPORT FidoExport
    FILE FidoTargets.cmake
    NAMESPACE Fido::
    DESTINATION ${CMAKE_INSTALL_PREFIX}/CMake
)

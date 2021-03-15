cmake_minimum_required(VERSION 3.15)

project(FIDO)

set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "")

include(GNUInstallDirs)

add_subdirectory(stdgpu)
add_subdirectory(ESSENS-master)
add_subdirectory(GRAAL)

install(EXPORT FidoExport
    FILE FidoTargets.cmake
    NAMESPACE Fido::
    DESTINATION ${CMAKE_INSTALL_PREFIX}/CMake
)

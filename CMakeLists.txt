cmake_minimum_required(VERSION 3.15)

project(FIDO)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

find_package(MPI REQUIRED)
find_package(Kokkos REQUIRED)
find_package(KokkosKernels REQUIRED)
find_package(resilience REQUIRED)
find_package(VeloC)
find_package(etl)

add_executable(fido src/main.cpp)
target_include_directories(fido PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(fido PRIVATE MPI::MPI_CXX)
target_link_libraries(fido PRIVATE Kokkos::kokkos)
target_link_libraries(fido PRIVATE Kokkos::kokkoskernels)
target_link_libraries(fido PRIVATE etl)
target_link_libraries(fido PRIVATE Kokkos::resilience)
target_link_libraries(fido PRIVATE veloc::veloc)

set_target_properties(fido PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
#set_target_properties(fido PROPERTIES
#    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
#)

install(TARGETS fido DESTINATION ${CMAKE_INSTALL_BINDIR})

#install(TARGETS fido
#    EXPORT fidoExport
#    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
#)
#install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/headers
#    EXPORT fidoExport
#    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
#    FILES_MATCHING PATTERN "*.hpp"
#)

#install(EXPORT fidoExport
#    FILE fidoTargets.cmake
#    NAMESPACE fido::
#    DESTINATION ${CMAKE_INSTALL_PREFIX}/CMake
#)

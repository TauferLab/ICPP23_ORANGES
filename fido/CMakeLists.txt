cmake_minimum_required(VERSION 3.15)

find_package(MPI REQUIRED)
find_package(Kokkos REQUIRED)
find_package(KokkosKernels REQUIRED)

set(FIDO_HEADERS
    headers/GDV_functions.hpp
    headers/bfs.hpp
    headers/class_definitions.hpp
    headers/combinations.hpp
    headers/print_disconnected_graph.hpp
)

add_executable(Fido parallel_graph_align/main_code_parallel.cpp ${FIDO_HEADERS})
target_include_directories(Fido PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/headers>
    $<INSTALL_INTERFACE:headers>
)
target_link_libraries(Fido PUBLIC ESSENS)
target_link_libraries(Fido PUBLIC MPI::MPI_CXX)
target_link_libraries(Fido PUBLIC Kokkos::kokkos)
target_link_libraries(Fido PUBLIC Kokkos::kokkoskernels)

set_target_properties(Fido PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
set_target_properties(Fido PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

install(TARGETS Fido
    EXPORT FidoExport
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/headers
    EXPORT FidoExport
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
)

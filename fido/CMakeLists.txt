cmake_minimum_required(VERSION 3.15)

find_package(MPI REQUIRED)
find_package(Kokkos REQUIRED)
find_package(KokkosKernels REQUIRED)
#find_package(resilience REQUIRED)
#find_package(VeloC)

set(TRACK_RUNTIME 1 CACHE STRING "Variable to define whether or not to record runtimes")
#set(THREAD_COUNT 4 CACHE STRING "Number of threads for shared memory kokkos paralellization")
set(SIM_MAT_CUT 100 CACHE STRING "Maximum similarity matrix size that will be printed to file")
set(USER_CHUNK_COUNT 64 CACHE STRING "Number of nodes assigned to a process at a time in a chunk")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/fido_compile_config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/include/fido_compile_config.h @ONLY)

set(FIDO_HEADERS
    include/GDV_functions.hpp
    include/bfs.hpp
    include/class_definitions.hpp
    include/combinations.hpp
    include/print_disconnected_graph.hpp
    include/fido_compile_config.h
)

add_executable(Fido src/main.cpp ${FIDO_HEADERS})
target_include_directories(Fido PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(Fido PRIVATE ESSENS)
target_link_libraries(Fido PRIVATE MPI::MPI_CXX)
target_link_libraries(Fido PRIVATE Kokkos::kokkos)
target_link_libraries(Fido PRIVATE Kokkos::kokkoskernels)
#target_link_libraries(Fido PRIVATE Kokkos::resilience)
#target_link_libraries(Fido PRIVATE veloc::veloc)

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
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include
    EXPORT FidoExport
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
)

cmake_minimum_required(VERSION 3.15)

find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_package(umpire REQUIRED)
find_package(RAJA REQUIRED)

set(FIDO_HEADERS
    headers/GDV_functions.hpp
    headers/GPUGDV_functions.hpp
    headers/bfs.hpp
    headers/gpubfs.hpp
    headers/class_definitions.hpp
    headers/combinations.hpp
    headers/combinations_raw.hpp
    headers/print_disconnected_graph.hpp
    headers/raw_vecs.hpp
)

set(CMAKE_CUDA_STANDARD 11)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda")

add_executable(fido main_code.cu ${FIDO_HEADERS})
target_include_directories(fido PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/headers>
    $<INSTALL_INTERFACE:headers>
)
target_link_libraries(fido PRIVATE Threads::Threads)
target_link_libraries(fido PUBLIC ESSENS)
target_link_libraries(fido PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(fido PUBLIC RAJA)
target_link_libraries(fido PUBLIC umpire)

set_target_properties(fido PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
set_target_properties(fido PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)
set_target_properties(fido PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

install(TARGETS fido
    EXPORT FidoExport
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/headers
    EXPORT FidoExport
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
)

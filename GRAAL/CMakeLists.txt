cmake_minimum_required(VERSION 3.15)

find_package(OpenMP REQUIRED)

set(FIDO_HEADERS
    headers/GDV_functions.hpp
    headers/bfs.hpp
    headers/class_definitions.hpp
    headers/combinations.hpp
    headers/print_disconnected_graph.hpp
)

add_executable(fido main_code.cpp ${FIDO_HEADERS})
target_include_directories(fido PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/headers>
    $<INSTALL_INTERFACE:headers>
)
target_link_libraries(fido PUBLIC ESSENS)
target_link_libraries(fido PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(fido PUBLIC stdgpu)
set_target_properties(fido PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
set_target_properties(fido PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

install(TARGETS fido
    EXPORT FidoExport
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/headers
    EXPORT FidoExport
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
)

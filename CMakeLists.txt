cmake_minimum_required(VERSION 3.15)

project(ORANGES)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")

find_package(OpenSSL REQUIRED)
find_package(MPI REQUIRED)
find_package(Kokkos REQUIRED)
find_package(KokkosKernels REQUIRED)
find_package(deduplicator REQUIRED)
#find_package(resilience REQUIRED)
#find_package(VeloC)
#find_package(etl)

file(GLOB SOURCES src/*.cpp)

add_executable(oranges ${SOURCES})

target_include_directories(oranges PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(oranges PRIVATE MPI::MPI_CXX)
target_link_libraries(oranges PRIVATE Kokkos::kokkos)
target_link_libraries(oranges PRIVATE Kokkos::kokkoskernels)
target_link_libraries(oranges PRIVATE Deduplicator::deduplicator)
target_link_libraries(oranges PRIVATE OpenSSL::SSL)
#target_link_libraries(oranges PRIVATE etl)
#target_link_libraries(oranges PRIVATE Kokkos::resilience)
#target_link_libraries(oranges PRIVATE veloc::veloc)

set_target_properties(oranges PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

install(TARGETS oranges DESTINATION ${CMAKE_INSTALL_BINDIR})


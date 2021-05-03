cmake_minimum_required(VERSION 3.15)

project(ESSENS)

# Change if versioning is added
set(ESSENS_VERSION "master")

set(ESSENS_SUBDIRS 
    Core/Basic_IO/PreProcess/Exceptions
    Core/Basic_IO/PreProcess/Level0
    Core/Basic_IO/PreProcess/Level1
    Core/Basic_IO/PreProcess/Level2
    Core/Basic_IO/Format/Level0
    Core/Basic_IO/Format/Input/Level1
    Core/Basic_IO/Format/Input/Level2
    Core/Basic_IO/Format/Output/Level1
    Core/Basic_IO/Format/Output/Level2
    Core/Basic_SetOps/Level0
    Core/Basic_SetOps/Level1
    Core/Basic_SetOps/Level2
    # Core/Basic_Traversal/Level0
    Core/Basic_Traversal/Level1
    Core/Basic_Traversal/Level2
    Basic_Analysis/Level0
)

add_library(ESSENS INTERFACE)
foreach(subdir ${ESSENS_SUBDIRS})
    target_include_directories(ESSENS INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${subdir}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
endforeach(subdir)

install(EXPORT EssensExport
    FILE EssensTargets.cmake
    NAMESPACE Essens::
    DESTINATION ${CMAKE_INSTALL_PREFIX}/CMake
)
install(TARGETS ESSENS
    EXPORT EssensExport
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
foreach(subdir ${ESSENS_SUBDIRS})
    message(STATUS "Subdir is ${subdir}")
    message(STATUS "Include Dir is ${CMAKE_INSTALL_INCLUDEDIR}")
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}
        EXPORT EssensExport
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING PATTERN "*.hpp"
    )
endforeach(subdir)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    EssensConfigVersion.cmake
    VERSION ${ESSENS_VERSION}
    COMPATIBILITY AnyNewerVersion
)
set(ESSENS_CMAKE_CONFIG_PATH "${CMAKE_INSTALL_PREFIX}/CMake")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/EssensConfig.cmake.in EssensConfig.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/EssensConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/EssensConfigVersion.cmake"
              DESTINATION ${CMAKE_INSTALL_PREFIX}/CMake
)

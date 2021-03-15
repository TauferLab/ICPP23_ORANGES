cmake_minimum_required(VERSION 3.15)

project(ESSENS)

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
        $<INSTALL_INTERFACE:${subdir}>
    )
endforeach(subdir)

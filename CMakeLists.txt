cmake_minimum_required(VERSION 2.8)

PROJECT(Endocardium_Points)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

add_executable(Endocardium_Points MACOSX_BUNDLE main.cpp )

target_link_libraries(Endocardium_Points ${VTK_LIBRARIES})

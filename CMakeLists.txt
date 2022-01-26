cmake_minimum_required(VERSION 3.14)

project(dsp VERSION 0.1 DESCRIPTION "DSP - A Modern C++ DSP library")

add_subdirectory(src ${CMAKE_SOURCE_DIR}/lib)
add_subdirectory(test ${CMAKE_SOURCE_DIR}/bin)
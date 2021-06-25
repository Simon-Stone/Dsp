cmake_minimum_required(VERSION 3.14)

project(dsp VERSION 0.1 DESCRIPTION "DSP - A Modern C++ DSP library")
list(APPEND CMAKE_MODULE_PATH "../scripts")

add_subdirectory(../src lib)
add_subdirectory(../test test)
enable_testing()

cmake_minimum_required(VERSION 3.14)

project(Dsp VERSION 0.1 DESCRIPTION "DSP - A Modern C++ DSP library")

add_subdirectory(src)
add_subdirectory(test)
enable_testing()
add_test(NAME example_test COMMAND UnitTests)
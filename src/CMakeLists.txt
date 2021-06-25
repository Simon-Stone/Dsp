# Some dependencies
find_package(FFTW3 REQUIRED)
find_package(FFTW3f REQUIRED)
find_package(FFTW3l REQUIRED)
include_directories(${FFTW3_INCLUDE_DIRS})
include_directories(${FFTW3f_INCLUDE_DIRS})
include_directories(${FFTW3l_INCLUDE_DIRS})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})


set(CMAKE_CXX_STANDARD 17)
add_library(dsp STATIC
        convert.cpp
        dsp.cpp
        fft.cpp
        filter.cpp
        Signal.cpp
        signals.cpp
        special.cpp
        utilities.cpp
        window.cpp
        )
target_include_directories(dsp PUBLIC "../include")
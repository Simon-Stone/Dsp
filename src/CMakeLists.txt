# Some dependencies
find_package(FFTW3)
if(FFTW3_FOUND)
    include_directories(${FFTW3_INCLUDE_DIRS})
endif()
find_package(FFTW3f)
if(FFTW3f_FOUND)
    include_directories(${FFTW3f_INCLUDE_DIRS})
endif()
find_package(FFTW3l)
if(FFTW3l_FOUND)
    include_directories(${FFTW3l_INCLUDE_DIRS})
endif()
find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR})
endif()

if(FFTW3_FOUND AND FFTW3f_FOUND AND FFTW3l_FOUND AND Boost_FOUND)
else()
add_definitions(-DZERO_DEPENDENCIES)
endif()

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
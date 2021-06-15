# Some dependencies
find_package(FFTW3 REQUIRED)
find_package(FFTW3f REQUIRED)
find_package(FFTW3l REQUIRED)
include_directories(${FFTW3_INCLUDE_DIR})
include_directories(${FFTW3f_INCLUDE_DIR})
include_directories(${FFTW3l_INCLUDE_DIR})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

set(CMAKE_CXX_STANDARD 17)
add_library(Dsp STATIC
        convert.cpp convert.h
        dsp.cpp dsp.h
        fft.cpp fft.h
        filter.cpp filter.h
        Signal.cpp Signal.h
        signals.cpp signals.h
        special.cpp special.h
        stats.h stats.h
        utilities.cpp utilities.h
        window.cpp window.h
        )

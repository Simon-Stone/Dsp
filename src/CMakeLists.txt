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
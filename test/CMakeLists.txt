# Some dependencies
find_package(FFTW3 REQUIRED)
find_package(FFTW3f REQUIRED)
find_package(FFTW3l REQUIRED)
include_directories(${FFTW3_INCLUDE_DIR})
include_directories(${FFTW3f_INCLUDE_DIR})
include_directories(${FFTW3l_INCLUDE_DIR})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Get Google Test
include(FetchContent)
FetchContent_Declare(
        googletest
        # Specify the commit you depend on and update it regularly.
        URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(CMAKE_CXX_STANDARD 17)
add_executable(UnitTests test.cpp)
target_link_directories(UnitTests PUBLIC ${FFTW3_LIBRARY_DIRS} ${FFTW3f_LIBRARY_DIRS} ${FFTW3l_LIBRARY_DIRS})
target_link_libraries(UnitTests PUBLIC gtest_main ${Boost_LIBRARIES} ${FFTW3_LIBRARIES} ${FFTW3f_LIBRARIES} ${FFTW3l_LIBRARIES} Dsp)

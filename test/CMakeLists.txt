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

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # GCC currently only supports the <execution> header through a third-party library
    find_package(TBB COMPONENTS tbbmalloc tbbmalloc_proxy tbb_preview)
endif()

set(CMAKE_CXX_STANDARD 17)

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

add_executable(UnitTests test.cpp)
list(APPEND libs "gtest_main" "dsp" "fftw3" "fftw3f" "fftw3l")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # GCC currently only supports the <execution> header through a third-party library
    list(APPEND libs "tbb")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Clang does not yet support the <execution> header
    message(FATAL_ERROR "Clang does not yet support the <execution> header. Please use GCC instead.")
endif()

target_link_directories(UnitTests PUBLIC ${FFTW3_LIBRARY_DIRS} ${FFTW3f_LIBRARY_DIRS} ${FFTW3l_LIBRARY_DIRS})
target_link_libraries(UnitTests PUBLIC ${libs})

include(GoogleTest)
gtest_discover_tests(UnitTests WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})



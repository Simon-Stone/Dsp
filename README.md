# Dsp
A digital signal processing library using Modern C++.

The library is supposed to be educational as much as it should be productive. It is therefore not highly optimized and there are definitely faster DSP libraries out there. But I am striving to keep the code clean, simple, and easy to read, maintain, and extend.

## Build instructions
### Dependencies
You can build the library without any external dependencies by defining a preprocessor constant ``ZERO_DEPENDENCIES``. You will miss out on some functionality and/or performance but at least you don't have to deal with missing libraries. 

If you choose to use the dependencies, you will need the following libraries:
1. FFTW
2. Boost::Math

**Important:** On Linux, you have to install Threading Building Blocks (TBB) because g++ relies on it for some of the C++17 parallel execution policies that are used.

#### Installing dependencies on Windows
I highly recommend using [vcpkg](https://vcpkg.info/) to install the required libraries. At long last, a package manager for C++ that works! 

TODO: Instructions using vcpkg

#### Installing dependencies on Linux
On Ubuntu:
``sudo apt-get install libboost-math-dev libfftw3-dev libtbb-dev``
Obviously, the package repository libraries will not be up-to-date, but if you are worried about that you probably know how to build them from soure.


### Windows using Visual Studio 2019
Open ``build/Dsp.sln``, choose your configuration, and build.


TODO: Instructions using CMake

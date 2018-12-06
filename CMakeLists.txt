cmake_minimum_required(VERSION 3.10)
project(stupid)

set(CMAKE_CXX_STANDARD 17)


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(lib/googletest)
include_directories(lib/googletest/googletest/include)

#include_directories(${GSL_INCLUDE_DIR})
#find_package(GSL REQUIRED)    # See below (2)


file(GLOB source src/*.cpp src/*.h)

add_library(stupid SHARED ${source})

target_link_libraries(stupid fftw3_omp fftw3 m gsl gslcblas)

target_include_directories(stupid PUBLIC include)


add_executable(test test/test.cpp)

target_link_libraries(test stupid gtest gtest_main pthread stdc++fs)

target_include_directories(test PUBLIC include src)
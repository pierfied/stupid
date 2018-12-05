cmake_minimum_required(VERSION 3.12)
project(stupid)

set(CMAKE_CXX_STANDARD 17)


file(GLOB source src/*.cpp src/*.h)

add_library(stupid SHARED ${source})

target_link_libraries(stupid m fftw3 gsl gslcblas)

target_include_directories(stupid PUBLIC include)


add_executable(test.exe test/test.cpp)

target_link_libraries(test.exe stupid gtest gtest_main pthread stdc++fs)

target_include_directories(test.exe PUBLIC include src)
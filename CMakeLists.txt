cmake_minimum_required(VERSION 3.8)
project(stupid)

set(CMAKE_CXX_STANDARD 17)


set(CMAKE_CXX_FLAGS_RELEASE "-O3")


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


file(GLOB source src/*.cpp src/*.h)

add_library(stupid SHARED ${source})

target_link_libraries(stupid m fftw3 fftw3_omp gsl gslcblas)

target_include_directories(stupid PUBLIC include)


add_executable(test.exe test/test.cpp)

target_link_libraries(test.exe stupid gtest gtest_main pthread stdc++fs)

target_include_directories(test.exe PUBLIC include src)
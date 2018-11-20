cmake_minimum_required(VERSION 3.12)
project(stupid)

set(CMAKE_CXX_STANDARD 17)

file(GLOB source src/*.h src/*.cpp)

add_library(stupid ${source})

add_executable(test ${source})
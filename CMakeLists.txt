cmake_minimum_required(VERSION 3.12)
project(stupid)

set(CMAKE_CXX_STANDARD 17)


file(GLOB source src/*.cpp src/*.h)

add_library(stupid ${source})

target_include_directories(stupid PUBLIC include)


file(GLOB test_source test/*.cpp test/*.h)

add_executable(test ${test_source})

target_link_libraries(test gtest stupid)

target_include_directories(test PUBLIC include src)
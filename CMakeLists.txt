cmake_minimum_required(VERSION 3.12)
project(stupid)

set(CMAKE_CXX_STANDARD 17)

add_library(stupid src/stupid.cpp src/stupid.h src/particle_list.cpp src/particle_list.h)

add_executable(test_particle_list src/particle_list.cpp src/particle_list.h)
cmake_minimum_required(VERSION 3.20)
project(untitled)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER /usr/bin/mpic++)

add_executable(1 main.cpp)
target_link_libraries(1 PUBLIC MPI::MPI_CXX)

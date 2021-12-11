cmake_minimum_required(VERSION 3.20)
project(lab2)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER /usr/bin/mpic++)

add_executable(lab2 main.cpp)
target_link_libraries(lab2 PUBLIC MPI::MPI_CXX)

cmake_minimum_required(VERSION 3.31)
project(LBM)

set(CMAKE_CXX_STANDARD 20)

include_directories(include)

add_executable(LBM main.cpp
        src/D2Q9.cpp
        src/D2Q9_cylinder.cpp
)

cmake_minimum_required(VERSION 3.31)
project(LBM)

set(CMAKE_CXX_STANDARD 20)

add_executable(LBM main.cpp
        src/D1Q3.cpp
        include/D1Q3.hpp
        src/D2Q9.cpp
        include/D2Q9.hpp
        src/D2Q9_cylinder.cpp
        include/D2Q9_cylinder.hpp)

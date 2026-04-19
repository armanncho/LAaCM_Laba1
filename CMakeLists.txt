cmake_minimum_required(VERSION 4.0)
project(linal1)

set(CMAKE_CXX_STANDARD 20)

add_executable(linal1 main.cpp
        matrix.cpp
        matrix.h
        gauss.cpp
        gauss.h
        lu.h
        lu.cpp
        utilits.cpp
        utilits.h)

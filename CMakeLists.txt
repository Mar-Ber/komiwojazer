cmake_minimum_required(VERSION 3.13)
project(komiwojazer)

set(CMAKE_CXX_STANDARD 14)

add_executable(komiwojazer main.cpp cmake-build-debug/ tsp.cpp tsp.hpp)
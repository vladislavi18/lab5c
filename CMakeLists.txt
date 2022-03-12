cmake_minimum_required(VERSION 3.20)
project(lab5c C)

set(CMAKE_C_STANDARD 11)

add_executable(lab5c main.c libs/algorithms/algorithms.c libs/algorithms/algorithms.h)

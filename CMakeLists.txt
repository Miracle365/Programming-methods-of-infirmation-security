cmake_minimum_required(VERSION 3.17)
project(Lab2 C)


link_libraries (crypto)


set(CMAKE_C_STANDARD 11)

add_executable(Lab2 main.c Mont_Curve.c Mont_Curve.h Point.c Point.h Parameters.h)
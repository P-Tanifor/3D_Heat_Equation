cmake_minimum_required(VERSION 3.16)
project(3D_Heat_Equation_)

set(CMAKE_CXX_STANDARD 14)

add_executable(3D_Heat_Equation_ mainThreadPool.cpp ThreadPool.h)
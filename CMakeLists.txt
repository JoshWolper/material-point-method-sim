cmake_minimum_required(VERSION 3.8)
project(MPM)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp main.cpp advection.cpp advection.h interpolation.h writeframe.h readfile.h mpmInitialize.h transfer.h global.h setBoundaryVelocity.h)
add_executable(MPM ${SOURCE_FILES})

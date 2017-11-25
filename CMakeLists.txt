cmake_minimum_required(VERSION 3.8)
project(MPM)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-O3 ")

set(SOURCE_FILES main.cpp main.cpp advection.h interpolation.h writeframe.h readfile.h mpmInitialize.h transfer.h global.h setBoundaryVelocity.h UpdateF.h SVD.h constitutiveModel.h Test/derivativeTest.cpp Test/derivativeTest.h Test/kernelTest.cpp Test/kernelTest.h interpolation.cpp constitutiveModel.cpp SVD.cpp computeMomentum.h AnalyticLevelSet.cpp AnalyticLevelSet.h AnalyticCollisionObject.cpp AnalyticCollisionObject.h)

add_executable(MPM ${SOURCE_FILES})

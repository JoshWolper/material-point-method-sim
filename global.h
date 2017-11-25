//
// Created by VElysianP on 10/1/2017.
//
#pragma once
#ifndef MATERIAL_POINT_METHOD_SIM_GLOBAL_H
#define MATERIAL_POINT_METHOD_SIM_GLOBAL_H

#include "Eigen/Eigen"

using namespace Eigen;

struct SVDResult
{
    Eigen::Matrix3f U;
    Eigen::Matrix3f SIGMA;
    Eigen::Matrix3f V;
};

struct SVDResultDouble
{
    Eigen::Matrix3d U;
    Eigen::Matrix3d SIGMA;
    Eigen::Matrix3d V;
};

struct Particle{
    float massP;
    float volumeP;
    Vector3f posP;
    Vector3f velP;
    Matrix3f BP;
    Matrix3f F;
    Matrix3f Fp;
    Matrix3f Fe;
};

struct GridAttr{
    float massG;
    Vector3f velG;//velocity at t_n+1
    Vector3f velGn;//velocity at t_n
    Vector3f force;
    Vector3f Xi;
};

struct GridInfo{
    float dx;
    int gridSize;
    int W;
    int H;
    int L;
};

#endif //MATERIAL_POINT_METHOD_SIM_GLOBAL_H

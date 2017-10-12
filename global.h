//
// Created by VElysianP on 10/1/2017.
//
#pragma once
#ifndef MATERIAL_POINT_METHOD_SIM_GLOBAL_H
#define MATERIAL_POINT_METHOD_SIM_GLOBAL_H

#include "Eigen/Eigen"

using namespace Eigen;

struct Particle{
    float massP;
    Vector3f posP;
    Vector3f velP;
    Matrix3f BP;
};

struct GridAttr{
    float massG;
    Vector3f velG;//velocity at t_n+1
    Vector3f velGn;//velocity at t_n
    Vector3f force;
};

struct GridInfo{
    float dx;
    int gridSize;
    int W;
    int H;
    int L;
};

#endif //MATERIAL_POINT_METHOD_SIM_GLOBAL_H

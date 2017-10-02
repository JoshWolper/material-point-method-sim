//
// Created by VElysianP on 10/1/2017.
//
#include "Eigen/Eigen"
#ifndef MATERIAL_POINT_METHOD_SIM_GLOBAL_H
#define MATERIAL_POINT_METHOD_SIM_GLOBAL_H

using namespace Eigen;

struct Particle{
    float massP;
    Vector3f posP;
    Vector3f velP;
    Matrix3f BP;
};

struct Grid{
    float massG;
    Vector3f velG;//velocity at t_n
    Vector3f velGn;//velocity at t_n+1
    Vector3f force;
};

#endif //MATERIAL_POINT_METHOD_SIM_GLOBAL_H

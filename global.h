//
// Created by VElysianP on 10/1/2017.
//

#ifndef MATERIAL_POINT_METHOD_SIM_GLOBAL_H
#define MATERIAL_POINT_METHOD_SIM_GLOBAL_H

#endif //MATERIAL_POINT_METHOD_SIM_GLOBAL_H

struct Particle{
    Vector3f posP;
    Vector3f velP;
    float massP;
    Matrix3f BP;
};

struct Grid{
    Vector3f posG;
    Vector3f velG;//vel n
    Vector3f velGN;//vel n+1
    float massG;
    Vector3f force;
};
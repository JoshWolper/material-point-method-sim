//
// Created by ziyinqu on 10/1/17.
//
#ifndef MPM_MPMINITIALIZE_H
#define MPM_MPMINITIALIZE_H

#include "readfile.h"
#include "global.h"
#include "Eigen/Eigen"

using namespace Eigen;

void mpmParticleInitialize(std::string filename, std::vector<Particle> &particles, float mass){
    std::vector<Vector3f> xp;
    readtxt(filename, xp);
    int part_number = xp.size();
    particles.resize(part_number);
    for (int i = 0; i < part_number; i++){
        particles[i].posP = xp[i];
        particles[i].velP = Vector3f::Zero();
        particles[i].massP = mass;
        particles[i].BP = Matrix3f::Zero();
    }
}

void mpmGridInitialize(std::vector<GridAttr> &gridAttr, GridInfo &gridInfo, Vector3i simArea, float dx){
    // Initialize grid information
    gridInfo.dx = dx;
    //TODO fix type?
    std::cout << simArea[0] << std::endl;
    std::cout << simArea(0) << std::endl;
    gridInfo.W = simArea[0]/dx;
    gridInfo.H = simArea[1]/dx;
    gridInfo.L = simArea[2]/dx;
    gridInfo.gridSize = gridInfo.W * gridInfo.H * gridInfo.L;
    // Initialize grid attribute
    gridAttr.resize(gridInfo.gridSize);
    for (int i = 0; i < gridInfo.gridSize; i++){
        gridAttr[i].massG = 0.f;
        gridAttr[i].force = Vector3f::Zero();
        gridAttr[i].velG = Vector3f::Zero();
        gridAttr[i].velGn = Vector3f::Zero();
    }

}

#endif //MPM_MPMINITIALIZE_H

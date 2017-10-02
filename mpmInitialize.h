//
// Created by ziyinqu on 10/1/17.
//
#include "readfile.h"
#include "global.h"
#include "Eigen/Eigen"
#ifndef MPM_MPMINITIALIZE_H
#define MPM_MPMINITIALIZE_H

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

void mpmGridInitialize(std::vector<Grid> &grid, float dx){
    int grid_size = 1/dx + 1;
    grid.resize(grid_size);
    for (int i = 0; i < grid_size; i++){
        grid[i].massG = 0.f;
        grid[i].force = Vector3f::Zero();
        grid[i].velG = Vector3f::Zero();
        grid[i].velGn = Vector3f::Zero();
    }

}

#endif //MPM_MPMINITIALIZE_H

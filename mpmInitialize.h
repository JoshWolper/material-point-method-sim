//
// Created by ziyinqu on 10/1/17.
//
#pragma once
#ifndef MPM_MPMINITIALIZE_H
#define MPM_MPMINITIALIZE_H

#include "readfile.h"
#include "global.h"
#include "Eigen/Eigen"

using namespace Eigen;

void mpmParticleInitialize(std::string filename, std::vector<Particle> &particles, float mass, float volume, Vector3f velocity){
    std::vector<Vector3f> xp;
    readtxt(filename, xp);
    int part_number = xp.size();
    std::cout << "INFO: >>>>>>>>>>>>>>> MPM Initialization <<<<<<<<<<<<<<<" << std::endl;
    std::cout << "INFO: Particle number is " << part_number << std::endl;
    particles.resize(part_number);
    for (int i = 0; i < part_number; i++){
        particles[i].posP = xp[i];
        particles[i].velP = velocity;
        particles[i].massP = mass;
        particles[i].volumeP = volume;
        particles[i].BP = Matrix3f::Zero();
        particles[i].F = Matrix3f::Identity();
        particles[i].Fe = Matrix3f::Identity();
        particles[i].Fp = Matrix3f::Identity();
    }
}

void mpmGridInitialize(std::vector<GridAttr> &gridAttr, GridInfo &gridInfo, Vector3i simArea, float dx){
    // Initialize grid information
    gridInfo.dx = dx;
    //TODO fix type?
    std::cout << "INFO: Simulation area is " << simArea[0] << "*" << simArea[1] << "*" << simArea[2]  << std::endl;
    gridInfo.W = simArea[0]/dx + 1;
    gridInfo.H = simArea[1]/dx + 1;
    gridInfo.L = simArea[2]/dx + 1;
    gridInfo.gridSize = gridInfo.W * gridInfo.H * gridInfo.L;
    std::cout << "INFO: Grid number is " << gridInfo.gridSize << std::endl;
    // Initialize grid attribute
    gridAttr.resize(gridInfo.gridSize);
    for (int i = 0; i < gridInfo.W; i++){
        for (int j = 0; j < gridInfo.H; j++){
            for (int k = 0; k < gridInfo.L; k++){
                int index = i * gridInfo.H * gridInfo.L + j * gridInfo.L + k;
                gridAttr[index].massG = 0;
                gridAttr[index].force = Vector3f::Zero();
                gridAttr[index].velG = Vector3f::Zero();
                gridAttr[index].velGn = Vector3f::Zero();
                gridAttr[index].Xi = Vector3f(i,j,k);
            }
        }
    }
    std::cout << "INFO: >>>>>>>>>>>>>>> MPM Initialization Ends <<<<<<<<<<<<<<< " << std::endl;
}

void mpmGridReinitialize(std::vector<GridAttr> &gridAttr, GridInfo &gridInfo){
    for (int i = 0; i < gridInfo.gridSize; i++){
        gridAttr[i].massG = 0;
        gridAttr[i].force = Vector3f::Zero();
        gridAttr[i].velG = Vector3f::Zero();
        gridAttr[i].velGn = Vector3f::Zero();
    }
}
#endif //MPM_MPMINITIALIZE_H

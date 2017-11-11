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

void mpmParticleInitialize(std::string filename, std::vector<Particle> &particles, float mass, float volume){
    std::vector<Vector3f> xp;
    readtxt(filename, xp);
    int part_number = xp.size();
    std::cout << "INFO: >>>>>>>>>>>>>>> MPM Initialization <<<<<<<<<<<<<<<" << std::endl;
    std::cout << "INFO: Particle number is " << part_number << std::endl;
    particles.resize(part_number);
    for (int i = 0; i < part_number; i++){
        particles[i].posP = xp[i];
        particles[i].velP = Vector3f::Zero();
        particles[i].massP = mass;
        particles[i].volumeP = volume;
        particles[i].BP = Matrix3f::Zero();
        particles[i].F = Matrix3f::Identity();
        particles[i].F = particles[i].F * 1.01;
        /*particles[i].F(0,0) = 1.0; //TODO: change this back to identity eventually!
        particles[i].F(0,1) = 2.0;
        particles[i].F(0,2) = 3.0;
        particles[i].F(1,0) = 4.0; //TODO: change this back to identity eventually!
        particles[i].F(1,1) = 5.0;
        particles[i].F(1,2) = 6.0;
        particles[i].F(2,0) = 7.0; //TODO: change this back to identity eventually!
        particles[i].F(2,1) = 8.0;
        particles[i].F(2,2) = 9.0;*/
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
    for (int i = 0; i < gridInfo.gridSize; i++){
        gridAttr[i].massG = 0.f;
        gridAttr[i].force = Vector3f::Zero();
        gridAttr[i].velG = Vector3f::Zero();
        gridAttr[i].velGn = Vector3f::Zero();
    }
    std::cout << "INFO: >>>>>>>>>>>>>>> MPM Initialization Ends <<<<<<<<<<<<<<< " << std::endl;
}

void mpmGridReinitialize(std::vector<GridAttr> &gridAttr, GridInfo &gridInfo){
    for (int i = 0; i < gridInfo.gridSize; i++){
        gridAttr[i].massG = 0.f;
        gridAttr[i].force = Vector3f::Zero();
        gridAttr[i].velG = Vector3f::Zero();
        gridAttr[i].velGn = Vector3f::Zero();
    }
}
#endif //MPM_MPMINITIALIZE_H

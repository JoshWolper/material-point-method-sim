//
// Created by ziyinqu on 10/1/17.
//
#pragma once
#ifndef MPM_ADVECTION_H
#define MPM_ADVECTION_H

#include "global.h"

using namespace std;

inline void addGravity(vector<GridAttr>& gridAttrs, vector<int> active_nodes, Vector3f gravity){
    for (int i = 0; i < active_nodes.size(); i++){
        int index = active_nodes[i];
        gridAttrs[index].force = gridAttrs[index].force + gridAttrs[index].massG * gravity ;
    }
}

inline void updateGridvelocity(vector<GridAttr>& gridAttrs, vector<int> active_nodes, float dt){
    for (int i = 0; i < active_nodes.size(); i++){
        int index = active_nodes[i];
        //Vector3f test = gridAttrs[index].force / gridAttrs[index].massG;
        //gridAttrs[index].velG = gridAttrs[index].velGn + dt * gridAttrs[index].force / gridAttrs[index].massG;
        gridAttrs[index].velG = Vector3f(0,-1,0);
    }
}

inline void addGridForces(vector<GridAttr>& gridAttrs, vector<Particle>& particles, GridInfo gridInfo, int energyDensityFunction) {

    //for each active grid node
    for (int i = 0; i < particles.size(); i++){
        //calculate force update
        float volume = particles[i].volumeP;
        Matrix3f defGrad = particles[i].F;
        float h = gridInfo.dx;

        //Calculate Piola Kirchoff Stress


        //Loop through each of the 27 grid nodes
            //Calculate gradW_ip for each of the grid nodes
            //Calculate f_i for each node
            //Add f_i to the node
    }
}

#endif //MPM_ADVECTION_H

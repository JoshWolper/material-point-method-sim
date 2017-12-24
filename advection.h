//
// Created by ziyinqu on 10/1/17.
//
#pragma once
#ifndef MPM_ADVECTION_H
#define MPM_ADVECTION_H

#include "global.h"
#include "constitutiveModel.h"

using namespace std;

void addGravity(vector<GridAttr>& gridAttrs, vector<int> active_nodes, Vector3f gravity){
    for (int i = 0; i < active_nodes.size(); i++){
        int index = active_nodes[i];
        gridAttrs[index].force = gridAttrs[index].force + gridAttrs[index].massG * gravity ;

    }
}

void addGridForces(vector<GridAttr>& gridAttrs, vector<Particle>& particles, GridInfo gridInfo, int energyDensityFunction) {
    //for each active grid node
    for (int i = 0; i < particles.size(); i++){
        //calculate force update
        float volume = particles[i].volumeP;
        Matrix3f defGrad = particles[i].F;
        Matrix3f Fp = particles[i].Fp;
        Matrix3f Fe = particles[i].Fe;
        float h = gridInfo.dx;

        //Calculate Piola Kirchoff Stress
        Matrix3f piola = Matrix3f::Ones();

        //TODO: small bottleneck
        switch(energyDensityFunction){
            case 0: corotatedPiola(defGrad, piola);
                break;
            case 1: neoHookeanPiola(defGrad, piola);
                break;
            case 2: stVernantPiola(defGrad, piola);
                break;
            case 3: snowPiola(defGrad, Fp, Fe, piola);
                break;
            default: corotatedPiola(defGrad, piola);
                break; //default should just be the corotated model
        }

        Vector3f pos = particles[i].posP / gridInfo.dx;
        Vector3i baseNode;
        Matrix3f wp = Matrix3f::Ones();
        Matrix3f dwp = Matrix3f::Ones();

        QuadraticInterpolation(pos, baseNode, wp, dwp); //get Wp and gradWp

        //Loop through each of the 27 grid nodes
        for(int r = 0; r < 3; r++){
            for(int s = 0; s < 3; s++){
                for(int t = 0; t < 3; t++) {

                    Vector3f gradWip;
                    float dwpx = dwp(0,r);
                    float wpy = wp(1,s);
                    float wpz = wp(2,t);
                    gradWip(0) = (1/h) * dwp(0,r) * wp(1, s) * wp(2, t); //calculate each component of gradWip
                    gradWip(1) = (1/h) * wp(0,r) * dwp(1, s) * wp(2, t);
                    gradWip(2) = (1/h) * wp(0,r) * wp(1, s) * dwp(2, t);

                    //cout << "Piola: " << piola << endl;
                    //cout << "defGrad: " << defGrad << endl;
                    //cout << "GradWip: " << gradWip << endl;

                    //TODO: bottleneck
                    Vector3f f_i = Vector3f::Zero();
                    if (energyDensityFunction != 3){
                        f_i = -1 * volume * piola * defGrad.transpose() * gradWip; //calc force update
                    }
                    else {
                        f_i = -1 * volume * piola * Fe.transpose() * gradWip; //calc force update
                    }


                    //cout << "f_i: " << f_i.transpose() << endl;

                    int x = baseNode(0) + r; //calculate the indeces of the node we're acting on
                    int y = baseNode(1) + s;
                    int z = baseNode(2) + t;
                    int index = (x * gridInfo.H * gridInfo.L) + (y * gridInfo.L) + z;
                    gridAttrs[index].force = gridAttrs[index].force + f_i; //add the force

                }
            }
        }

    }
}

void updateGridvelocity(vector<GridAttr>& gridAttrs, vector<int> active_nodes, float dt){
    for (int i = 0; i < active_nodes.size(); i++){
        int index = active_nodes[i];
        Vector3f deltav = dt * gridAttrs[index].force / gridAttrs[index].massG;
        gridAttrs[index].velG = gridAttrs[index].velGn + deltav;
    }
}


#endif //MPM_ADVECTION_H

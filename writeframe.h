//
// Created by ziyinqu on 10/1/17.
//
#pragma once
#ifndef MPM_WRITEFRAME_H
#define MPM_WRITEFRAME_H

#include <string>
#include "global.h"

using namespace std;
void saveStep(vector<Particle> particles, int step){

    //Write to an object file
    ofstream outfile;

    //later we'll want a second parameter, timeStep to make it spit out differently named files!
    string filename = "Output/step" + to_string(step) + ".obj";

    outfile.open(filename);

    for(int i = 0; i < particles.size(); i++){
        outfile << "v " << particles[i].posP[0] << " " << particles[i].posP[1] << " " << particles[i].posP[2] << "\n";
    }
    outfile.close();

    return;
}

void saveFrame(vector<Particle> particles, int frame){

    //Write to an object file
    ofstream outfile;

    //later we'll want a second parameter, timeStep to make it spit out differently named files!
    string filename = "Output/frame" + to_string(frame) + ".obj";

    outfile.open(filename);

    for(int i = 0; i < particles.size(); i++){
        outfile << "v " << particles[i].posP[0] << " " << particles[i].posP[1] << " " << particles[i].posP[2] << "\n";
    }
    outfile.close();

    return;
}

#endif //MPM_WRITEFRAME_H

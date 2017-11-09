//
// Created by ziyinqu on 10/1/17.
//
#pragma once
#ifndef MPM_WRITEFRAME_H
#define MPM_WRITEFRAME_H

#include <string>

using namespace std;
void saveStep(vector<Vector3f> points, int step){

    //Write to an object file
    ofstream outfile;

    //later we'll want a second parameter, timeStep to make it spit out differently named files!
    string filename = "Output/step" + to_string(step) + ".obj";

    outfile.open(filename);

    for(int i = 0; i < points.size(); i++){
        double point [3];
        point[0] = points[i][0];
        point[1] = points[i][1];
        point[2] = points[i][2];
        outfile << "v " << point[0] << " " << point[1] << " " << point[2] << "\n";
    }
    outfile.close();

    return;
}

void saveFrame(vector<Vector3f> points, int frame){

    //Write to an object file
    ofstream outfile;

    //later we'll want a second parameter, timeStep to make it spit out differently named files!
    string filename = "Output/frame" + to_string(frame) + ".obj";

    outfile.open(filename);

    for(int i = 0; i < points.size(); i++){
        double point [3];
        point[0] = points[i][0];
        point[1] = points[i][1];
        point[2] = points[i][2];
        outfile << "v " << point[0] << " " << point[1] << " " << point[2] << "\n";
    }
    outfile.close();

    return;
}

#endif //MPM_WRITEFRAME_H

//
// Created by ziyinqu on 10/1/17.
//
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Eigen/Eigen"

#ifndef MPM_READFILE_H
#define MPM_READFILE_H

using namespace Eigen;

void readtxt(std::string filename, std::vector<Vector3f>& xp){
    std::ifstream inputfile;
    inputfile.open(filename);
    if(!inputfile){
        std::cerr << "Unable to open " << filename << "!" << std::endl;
        exit(1);
    }
    std::string line;
    while(std::getline(inputfile, line)){
        std::stringstream ss(line);
        Vector3f thisXp;
        for (int i = 0; i < 3; i++){
            ss >> thisXp(i);
            //TODO assert data
            //assert(thisXp(i) > 0 && thisXp(i) < 1);
        }
        xp.push_back(thisXp);
    }
}
#endif //MPM_READFILE_H

//
// Created by ziyinqu on 10/1/17.
//

#include "readfile.h"

int main(){
    std::string file = "sparseDragonSamples.txt";
    std::vector<Vector3f> xp;
    std::vector<Vector3f> vp;
    readtxt(file, xp);
    int size = xp.size();
    vp.resize(size);
    return 0;
}
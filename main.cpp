//
// Created by ziyinqu on 10/1/17.
//
#include <string>
#include "mpmInitialize.h"
#include "global.h"

int main(){

    // particles attributes initialize
    std::string filename = "sparseDragonSamples.txt";
    std::vector<Particle> particles;
    float mass = 1.f;
    mpmParticleInitialize(filename, particles, mass);

    // grid attributes initialize
    float dx = 0.02f;
    std::vector<Grid> grid;
    mpmGridInitialize(grid, dx);

    return 0;
}
//
// Created by ziyinqu on 10/1/17.
//
#include <string>
#include "mpmInitialize.h"
#include "transfer.h"

int main(){

    // MPM simulation parameters setting up
    float dt = 0.02f;
    float alpha = 0.95;

    // particles attributes initialize
    float mass = 0.1f;
    std::string filename = "sparseDragonSamples.txt";
    std::vector<Particle> particles;
    mpmParticleInitialize(filename, particles, mass);

    // grid attributes initialize
    float dx = 0.02f;
    Vector3i simArea = Vector3i::Ones();
    std::vector<GridAttr> gridAttrs;
    GridInfo gridInfo;
    mpmGridInitialize(gridAttrs, gridInfo, simArea, dx);

    // transfer from Particles to Grid
    transferP2G(particles, gridAttrs, gridInfo);

    //TODO add advection part here

    //TODO add boudnary collision here

    //TODO update deformation gradient here

    // transfer from Grid to particles
    transferG2P(particles, gridAttrs, gridInfo, dt, alpha);
    return 0;
}
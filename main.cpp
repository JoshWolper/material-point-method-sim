//
// Created by ziyinqu on 10/1/17.
//
#include <string>
#include "mpmInitialize.h"
#include "transfer.h"
#include "advection.h"
#include "writeframe.h"

int main(){

    // MPM simulation parameters setting up
    float dt = 0.02f;
    float alpha = 0.95;
    Vector3f gravity = Vector3f(0, 0, -9.8f);
    std::vector<int> active_nodes;

    // particles attributes initialize
    float mass = 0.1f;
    std::string filename = "Models/sparseDragonSamples.txt";
    std::vector<Particle> particles;
    mpmParticleInitialize(filename, particles, mass);

    // grid attributes initialize
    float dx = 0.02f;
    Vector3i simArea = Vector3i::Ones();
    std::vector<GridAttr> gridAttrs;
    GridInfo gridInfo;
    mpmGridInitialize(gridAttrs, gridInfo, simArea, dx);

    // transfer from Particles to Grid
    transferP2G(particles, gridAttrs, gridInfo, active_nodes);

    // advection part, add forces and update grid velocity
    addGravity(gridAttrs, active_nodes, gravity);

    updateGridvelocity(gridAttrs, active_nodes, dt);

    //TODO add boudnary collision here

    //TODO update deformation gradient here

    // transfer from Grid to particles
    transferG2P(particles, gridAttrs, gridInfo, dt, alpha);

    //Save the particles!
    vector<Vector3f> points;
    for(int i = 0; i < particles.size(); i++){
	    points.push_back(particles[i].posP); //add the positions to the points list!
    }
    saveParticles(points); //call the saving routine

    return 0;
}

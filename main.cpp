//
// Created by ziyinqu on 10/1/17.
//
#include <string>
#include "mpmInitialize.h"
#include "transfer.h"

void saveParticles(vector<Vector3f> points);

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

    //Save the particles!
    vector<Vector3f> points;
    for(int i = 0; i < particles.size(); i++){
	points.push_back(particles[i].posP); //add the positions to the points list!
    }
    saveParticles(points); //call the saving routine

    return 0;
}

void saveParticles(vector<Vector3f> points){

    //Write to an object file
    ofstream outfile;
    
    //later we'll want a second parameter, timeStep to make it spit out diferently named files!
    //String filename = "sim_t=" + itos(

    outfile.open("simOut.obj");

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

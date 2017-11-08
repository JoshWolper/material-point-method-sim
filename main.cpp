//
// Created by ziyinqu on 10/1/17.
//
#include <string>
#include "mpmInitialize.h"
#include "transfer.h"
#include "advection.h"
#include "writeframe.h"
#include "setBoundaryVelocity.h"
#include "UpdateF.h"

int main(){

    // MPM simulation parameters setting up
    float dt = 0.02f;
    float alpha = 0;
    Vector3f gravity = Vector3f(0, -9.8, 0);

    // particles attributes initialize
    float mass = 0.1f;
    float volume = 0.1f;
    std::string filename = "Models/sparseDragonTranslated.obj";
    std::vector<Particle> particles;
    mpmParticleInitialize(filename, particles, mass, volume);

    // grid attributes initialize
    float dx = 0.04f;
    Vector3i simArea = Vector3i::Ones();
    std::vector<GridAttr> gridAttrs;
    GridInfo gridInfo;
    mpmGridInitialize(gridAttrs, gridInfo, simArea, dx);

    // start simulation
    int step = 0;
    cout << "INFO: >>>>>>>>>>>>>>> Simulation Start! <<<<<<<<<<<<<<< " << endl;
    while (step != 15) {
        cout << "INFO: Current simulation step is " << step << endl;
        mpmGridReinitialize(gridAttrs, gridInfo);
        std::vector<int> active_nodes;
        // transfer from Particles to Grid
        transferP2G(particles, gridAttrs, gridInfo, active_nodes);

        // advection part, add forces and update grid velocity
        addGravity(gridAttrs, active_nodes, gravity);

        //Add external forces based on our defined energy density function
        int energyDensityFunction = 0; //define which density function we wish to use!
        addGridForces(gridAttrs, particles, gridInfo, energyDensityFunction);

        updateGridvelocity(gridAttrs, active_nodes, dt);

        // boundary collision
        setBoundaryVelocity(gridAttrs, gridInfo);

        //update deformation gradient here
        //TODO add loop over particles in UpdateF
        UpdateF(dt, gridInfo, gridAttrs, particles[0]);

        // transfer from Grid to particles
        transferG2P(particles, gridAttrs, gridInfo, dt, alpha);

        //Save the particles!
        vector<Vector3f> points;
        for (int i = 0; i < particles.size(); i++) {
            points.push_back(particles[i].posP); //add the positions to the points list!
        }
        saveParticles(points, step); //call the saving routine
        step = step + 1;
    }
    return 0;
}

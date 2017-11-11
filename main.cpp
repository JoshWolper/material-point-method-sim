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
#include "Test/kernelTest.h"
#include <ctime>

#define SANITYCHECK false

int main(){

    if (SANITYCHECK){
        quadraticTest();
    }
    else {
        // MPM simulation parameters setting up
        float dt = 1e-3f; //50 FPS
        float frameRate = 24;
        int stepsPerFrame = (int)ceil(1 / (dt / (1 / frameRate))); //calculate how many steps per frame we should have based on our desired frame rate and dt!
        float alpha = 0;
        Vector3f gravity = Vector3f(0, 0, 0);

        // particles attributes initialize
        float density = 1.0f;
        int numPoints = 3653; //points in sparse cube
        float mass = 10;
        float volume = mass/density;
        //std::string filename = "../Models/VeryDenseCube.obj";
        std::string filename = "Models/newSparseCube_Nov9.obj";
        std::vector<Particle> particles;
        mpmParticleInitialize(filename, particles, mass, volume);

        // grid attributes initialize
        float dx = 0.02f;
        Vector3i simArea = Vector3i::Ones();
        std::vector<GridAttr> gridAttrs;
        GridInfo gridInfo;
        mpmGridInitialize(gridAttrs, gridInfo, simArea, dx);

        // start simulation
        int step = 0;
        int frame = 0;
        cout << "INFO: >>>>>>>>>>>>>>> Simulation Start! <<<<<<<<<<<<<<< " << endl;
        while (step != 1000) {
            // set timer
            std::clock_t start, p2gstart, g2pstart, updatestart, forcestart, updatefstart;
            double duration, p2gduration, g2pduration, updateduration, forceduration, updatefduration;
            start = std::clock();

            cout << "INFO: Current simulation step is " << step << endl;
            mpmGridReinitialize(gridAttrs, gridInfo);
            std::vector<int> active_nodes;
            // transfer from Particles to Grid
            p2gstart = std::clock();
            transferP2G(particles, gridAttrs, gridInfo, active_nodes);
            p2gduration = ( std::clock() - p2gstart ) / (double) CLOCKS_PER_SEC;

            // advection part, add forces and update grid velocity
            //addGravity(gridAttrs, active_nodes, gravity);

            //Add external forces based on our defined energy density function
            forcestart = std::clock();
            int energyDensityFunction = 0; //define which density function we wish to use!
            addGridForces(gridAttrs, particles, gridInfo, energyDensityFunction);
            forceduration = ( std::clock() - forcestart ) / (double) CLOCKS_PER_SEC;

            updatestart = std::clock();
            updateGridvelocity(gridAttrs, active_nodes, dt);
            updateduration = ( std::clock() - updatestart ) / (double) CLOCKS_PER_SEC;

            // boundary collision
            setBoundaryVelocity(gridAttrs, gridInfo);

            //update deformation gradient here
            updatefstart = std::clock();
            UpdateF(dt, gridInfo, gridAttrs, particles);
            updatefduration = ( std::clock() - updatefstart ) / (double) CLOCKS_PER_SEC;

            // transfer from Grid to particles
            g2pstart = std::clock();
            transferG2P(particles, gridAttrs, gridInfo, dt, alpha);
            g2pduration = ( std::clock() - g2pstart ) / (double) CLOCKS_PER_SEC;


            if (step % stepsPerFrame == 0) {
                //Save the particles!
                cout << "INFO: Current Frame is " << frame << endl;
                saveFrame(particles, frame); //call the saving routine
                frame = frame + 1;
            }
                //Save the particles!
//                saveStep(particles, step); //call the saving routine

            step = step + 1;
            // output time calculation
            duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
            cout << "      P2G time is: " << p2gduration << endl;
            cout << "      UPV time is: " << updateduration << endl;
            cout << "      FRC time is: " << forceduration << endl;
            cout << "      UPF time is: " << updatefduration << endl;
            cout << "      G2P time is: " << g2pduration << endl;
            cout << "      overall time is: " << duration << endl;
        }
    }

    return 0;
}

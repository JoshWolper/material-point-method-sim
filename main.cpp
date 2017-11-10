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

#define SANITYCHECK true
int main(){

    if (SANITYCHECK){
        quadraticTest();
    }
    else {
        // MPM simulation parameters setting up
        float dt = 1e-3f; //50 FPS
        float frameRate = 24;
        int stepsPerFrame = ceil(1 / (dt / (1 /
                                            frameRate))); //calculate how many steps per frame we should have based on our desired frame rate and dt!
        float alpha = 0;
        Vector3f gravity = Vector3f(0, 0, 0);

        // particles attributes initialize
        float density = 1.0f;
        int numPoints = 3653; //points in sparse cube
        float mass = (((float) 1 / (float) 27) * density) / numPoints;
        float volume = ((float) 1 / (float) 27) / numPoints;
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
            UpdateF(dt, gridInfo, gridAttrs, particles);

            // transfer from Grid to particles
            transferG2P(particles, gridAttrs, gridInfo, dt, alpha);

            if (step % stepsPerFrame == 0) {
                //Save the particles!
                vector<Vector3f> points;
                for (int i = 0; i < particles.size(); i++) {
                    points.push_back(particles[i].posP); //add the positions to the points list!
                }
                saveFrame(points, frame); //call the saving routine
                frame = frame + 1;
            }

            step = step + 1;
        }
    }

    return 0;
}

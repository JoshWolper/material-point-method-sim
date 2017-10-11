//
// Created by ziyinqu on 10/1/17.
//

#ifndef MPM_WRITEFRAME_H
#define MPM_WRITEFRAME_H

void saveParticles(vector<Vector3f> points){

    //Write to an object file
    ofstream outfile;

    //later we'll want a second parameter, timeStep to make it spit out differently named files!
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

#endif //MPM_WRITEFRAME_H

//
// Created by Josh Wolper on 10/11/17.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <ctime>
#include <Eigen/Eigen>
#include </../../../../../../usr/include/eigen3/Eigen/Eigen>

main(int argc, char* argv[]) {

    vector<Vector3f> points;

    if (argc != 2) { // argc should be 2 for correct execution
        // We print argv[0] assuming it is the program name
        cout << "usage: " << argv[0] << " <targetModel.txt>\n";
        return 0;
    } else {
        // We assume argv[1] is a filename to open
        FILE *afile = fopen(argv[1], "r");
        // Always check to see if file opening succeeded
        if (afile == NULL) {
            cout << "Could not open first file\n";
            return 0;
        }
    }

    while ( 1 ){
        int ret = fscanf(afile, "%f %f %f\n", &x, &y, &z);
        if(ret == EOF){
            break;
        }
        else{
            points.push_back(Vector3f(x, y, z));
        }
    }


    Vector3f min = Vector3f(-100,-100,-100);
    Vector3f max = Vector3f(100,100,100);

    //Now find the min and max of the model
    for(int i = 0; i < points.size(); i++){

        Vector3f currPoint = points[i];

        if((min[0] < currPoint[0]) && (min[1] < currPoint[1]) && (min[2] < currPoint[2])){

            min = currPoint;

        }
        else if( (max[0] > currPoint[0]) && (max[1] > currPoint[1]) && (max[2] > currPoint[2])){

            max = currPoint;

        }

    }

    Vector3f modelDimensions = max - min;

    float maxDim = -1;
    for(int i = 0; i < 3; i++){

        if(modelDimensions[i] > maxDim){

            maxDim = modelDimensions[i];

        }

    }


    float dx = 0.5;
    float dy = 0.6;

    for(int i = 0; i < points.size(); i++){

        points[i] = points[i] / maxDim; //scale all points down by the maximum dimension, should make every dimension now between 0 and 1!

        points[i] = Vector3f(points[i][0] + dx, points[i][1] + dy, points[i][2]); //add the translation based on dy and dx!

    }

    saveSamples(points);

    return;

}

void saveSamples(vector<Vector3f> points){

    //Write to a binary file and a text file
    ofstream outfile;

    outfile.open("translatedAndScaledModel.txt");

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





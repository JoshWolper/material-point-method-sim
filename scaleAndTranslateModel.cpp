//
// Created by Josh Wolper on 10/11/17.
//

#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointCount.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <ctime>
#include </../../../../../../usr/include/eigen3/Eigen/Eigen>

using namespace std;
using namespace openvdb;

void saveSamples(vector<Vec3d> points);

main(int argc, char* argv[]) {

    vector<Vec3d> points;

    if (argc != 2) { // argc should be 2 for correct execution
        // We print argv[0] assuming it is the program name
        cout << "usage: " << argv[0] << " <targetModel.txt>\n";
        return 0;
    } else {
        // We assume argv[1] is a filename to open
        FILE * afile = fopen(argv[1], "r");
        // Always check to see if file opening succeeded
        if (afile == NULL) {
            cout << "Could not open first file\n";
            return 0;
        }

	float x,y,z;

    	while ( 1 ){
        	int ret = fscanf(afile, "%f %f %f\n", &x, &y, &z);
        	if(ret == EOF){
            		break;
        	}
        	else{
            		points.push_back(Vec3d(x, y, z));
        	}
    	}


    	Vec3d min = Vec3d(-100,-100,-100);
    	Vec3d max = Vec3d(100,100,100);

    	//Now find the min and max of the model
    	for(int i = 0; i < points.size(); i++){

        	Vec3d currPoint = points[i];

       	 	if((min.x() > currPoint.x()) && (min.y() > currPoint.y()) && (min.z() > currPoint.z())){

            		min = currPoint;

        	}
        	else if( (max.x() < currPoint.x()) && (max.y() < currPoint.y()) && (max.z() < currPoint.z())){

            		max = currPoint;

        	}

    	}

    	Vec3d modelDimensions = max - min;

	double yTranslate = modelDimensions.y()/2;

	//Translate to 0,0,0
	for(int i = 0; i < points.size(); i++){

        	//points[i] = Vec3d(points[i].x(), points[i].y() - yTranslate, points[i].z()); //shift down in YDimension
    	}
	
    	double maxDim = -1;
    	for(int i = 0; i < 3; i++){

        	if(modelDimensions[i] > maxDim){

            		maxDim = modelDimensions[i]; //only scale by a proportion of this

        	}

    	}

	//float scaleFactor = maxDim / 5;
	float scaleFactor = 22;

    	double dx = 0.5;
    	double dy = 0.2;
	double dz = 0.5;

    	for(int i = 0; i < points.size(); i++){

        	points[i] = points[i] / scaleFactor; //scale all points down by the maximum dimension, should make every dimension now between 0 and 1!

        	points[i] = Vec3d(points[i].x() + dx, points[i].y() + dy, points[i].z() + dz); //add the translation based on dy and dx!

    	}

    	saveSamples(points);

    }

    

    return 0;

}

void saveSamples(vector<Vec3d> points){

    //Write to a binary file and a text file
    ofstream outfile;

    outfile.open("translatedAndScaledModel.txt");

    for(int i = 0; i < points.size(); i++){
        double point [3];
        point[0] = points[i].x();
        point[1] = points[i].y();
        point[2] = points[i].z();
        outfile << "v " << point[0] << " " << point[1] << " " << point[2] << "\n";
    }
    outfile.close();

    return;
}





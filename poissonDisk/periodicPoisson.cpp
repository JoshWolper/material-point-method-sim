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
#include </../../../../../../usr/include/eigen3/Eigen/Dense>

using namespace std;
using namespace openvdb;

void saveSamples(vector<Vec3d> samples);

main(int argc, char* argv[]){
	
	openvdb::initialize(); //init at beginning so it knows about all our vars
	
	FloatGrid::Ptr grid;

	vector<Vec3d> sampleCube;

	vector<Vec3d> samples;

	if ( argc != 5 ){ // argc should be 2 for correct execution
    		// We print argv[0] assuming it is the program name
    		cout<<"usage: "<< argv[0] <<" <levelSetFile.vdb> <densePointsFile.txt> <rValueDesired> <rValueSampleCube";
		return 0;
  	}
	else {
    		// We assume argv[1] is a filename to open
    		FILE * afile = fopen( argv[1], "r" );
    		// Always check to see if file opening succeeded
    		if ( afile == NULL ){
      			cout<< "Could not open first file\n";
			return 0;
    		}
		FILE * afile2 = fopen( argv[2], "r");
		if(afile2 == NULL){
			cout << "Could not open second file\n";
			return 0;
		}
		//Here we'll load up the VDB level set
		io::File file(argv[1]);
		file.open();
		GridBase::Ptr baseGrid;
		for(io::File::NameIterator nameIter = file.beginName(); nameIter != file.endName(); ++nameIter){
			baseGrid = file.readGrid(nameIter.gridName());
		}
		file.close();
		grid = gridPtrCast<FloatGrid>(baseGrid);

		//Initialize a grid sampler so we can find the values of continuous points sampled in the grid!		
		tools::GridSampler<FloatTree, tools::BoxSampler> interpolator(grid->constTree(), grid->transform());

		float rDesired = stof(argv[3]);
		float rCube = stof(argv[4]);

		float scaleFactor = rDesired / rCube;

		float x,y,z;

		while ( 1 ){
        		int ret = fscanf(afile2, "%f %f %f\n", &x, &y, &z);
			if(ret == EOF){
				break;
			}
			else{
				sampleCube.push_back(Vec3d(x * scaleFactor, y * scaleFactor, z * scaleFactor));
			}
		}

		CoordBBox boundbox;
		Coord min, max;
		boundbox = grid->evalActiveVoxelBoundingBox();
		min = boundbox.min();
		max = boundbox.max();

		cout << "Model Min: " << min.x() << " " << min.y() << " " << min.z() << "\n";
		cout << "Model Max: " << max.x() << " " << max.y() << " " << max.z() << "\n";

		cout << "Finished loading everything up! Beginning periodic tiling!\n";
	
		float cubeDim = scaleFactor; //1x1x1 originally so it is just the scalefactor
		Vec3d cubeMin = Vec3d( -0.5 * scaleFactor, -0.5 * scaleFactor, -0.5 * scaleFactor);
		Vec3d cubeMax = Vec3d( 0.5 * scaleFactor, 0.5 * scaleFactor, 0.5 * scaleFactor);

		int xDim, yDim, zDim;
		xDim = (int)ceil(abs(max.x() - min.x()) / cubeDim); //figure out how many cubes we need in each dimension!
		yDim = (int)ceil(abs(max.y() - min.y()) / cubeDim);
		zDim = (int)ceil(abs(max.z() - min.z()) / cubeDim);

		Vec3d cubeCenter = Vec3d(min.x() + (0.5 * scaleFactor), min.y() + (0.5 * scaleFactor), min.z() + (0.5 * scaleFactor));
		cout << "Initial cube center: " << cubeCenter.x() << " "<< cubeCenter.y() << " " << cubeCenter.z() << "\n";

		for(int i = 0; i < yDim; i++){
			for(int j = 0; j < xDim; j++){
				for(int k = 0; k < zDim; k++){
					cout << "Current indeces: " << i << " " << j << " " << k << "\n";
					//Get new cube center! (we're tiling the cube through the model bounding box)
					cubeCenter = Vec3d(min.x() + (j * scaleFactor), min.y() + (i * scaleFactor), min.z() + (k * scaleFactor));
					
					cout << "Current cube center: " << cubeCenter.x() << " "<< cubeCenter.y() << " " << cubeCenter.z() << "\n";

					//Now for each point we have to add the center and then check whether that point is inside the mesh, if so keep it!
					for(int p = 0; p < sampleCube.size(); p++){
						Vec3d newPoint = sampleCube[p] + cubeCenter;
						if(interpolator.wsSample(newPoint) < 0){
							samples.push_back(newPoint); //inside mesh! keep this sample
						}
					}
				}
			}
		}
		
		cout << "Number of samples: " << samples.size() << "\n";
		saveSamples(samples);
	}

	return 0;

}

void saveSamples(vector<Vec3d> samples){
	
	//Write to a binary file and a text file
	ofstream outfile;
	ofstream outfile2;
	outfile.open("novelSamples.dat", ios::binary | ios::out);
	outfile2.open("novelSamples.txt");
	for(int i = 0; i < samples.size(); i++){
		double point [3];
		point[0] = samples[i].x();
		point[1] = samples[i].y();
		point[2] = samples[i].z();
		//outfile.write(point, 3*sizeof(double)); 
		outfile2 << point[0] << " " << point[1] << " " << point[2] << "\n";
	}
	outfile.close();
	outfile2.close();
	
	//Also write to a VDB file
	openvdb::initialize();

	vector<Vec3R> positions;
	Vec3d currSample;
	while(!samples.empty()){
		currSample = samples.back(); //get the last sample in the list
		positions.push_back(Vec3R(currSample[0], currSample[1], currSample[2])); //add to the positions list
		samples.pop_back(); //pop the final sample in the list since we've now added it!
	}

	//Set up a point grid and write to a file so we can view the output points!
	openvdb::points::PointAttributeVector<Vec3R> positionsWrapper(positions);	
	int pointsPerVoxel = 8;
	float voxelSize = points::computeVoxelSize(positionsWrapper, pointsPerVoxel);
	openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(voxelSize);
	openvdb::points::PointDataGrid::Ptr grid = points::createPointDataGrid<points::NullCodec, points::PointDataGrid>(positions, *transform);
	grid->setName("Points");
	io::File file("denseModel.vdb");
	GridPtrVec grids;
	grids.push_back(grid);
	file.write(grids);
	file.close();	


	return;
}


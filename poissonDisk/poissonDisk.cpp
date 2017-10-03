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

vector<float> randPointInBBox(Vec3d min, Vec3d max);
vector<int> backgroundGridIndeces(Vec3d point, Vec3d min, float bucketSize);
Vec3d getPointInAnnulus(Vec3d center, float r);
int checkNeighborCells(vector<int> currIndex, Vec3d currPoint, vector<vector<vector<int>>> backGrid, vector<Vec3d> samples, float r);
int distOk(Vec3d p1, Vec3d p2, float r);
void saveSamples(vector<Vec3d> samples);

main(int argc, char* argv[]){
	
	srand(static_cast <unsigned> (time(0))); //seed the rng	

	openvdb::initialize(); //init at beginning so it knows about all our vars
	
	FloatGrid::Ptr grid;

	if ( argc != 5 ){ // argc should be 2 for correct execution
    		// We print argv[0] assuming it is the program name
    		cout<<"usage: "<< argv[0] <<" <filename> <scaleFactor> <rValue> <kValue>";
  	}
	else {
    		// We assume argv[1] is a filename to open
    		FILE * afile = fopen( argv[1], "r" );
    		// Always check to see if file opening succeeded
    		if ( afile == NULL ){
      			cout<< "Could not open file\n";
			return 0;
    		}
		else {
			//Here we'll load up the VDB level set
			io::File file(argv[1]);
			file.open();
			GridBase::Ptr baseGrid;
			for(io::File::NameIterator nameIter = file.beginName(); nameIter != file.endName(); ++nameIter){
				baseGrid = file.readGrid(nameIter.gridName());
			}
			file.close();
			grid = gridPtrCast<FloatGrid>(baseGrid);
   		}

		cout << "Finished opening level set.\n";
		cout << "Poisson Disk Algorithm begin.\n";
		
		//Poisson Disk Algorithm Begin!

		int n = 3; //define our dimension (3D)
		float r = stof(argv[3]); //define minimum distance between samples
		int kParam = stoi(argv[4]); //define number of samples to take before rejecting current sample

		//Initialize 3D array for the background grid for storing samples, store all -1 at first
		//Cell size is bounded by r/sqrt(n) sp if we get the bounding box for the level set we can determine how big the array must be in each dimension

		//Get the bounding box and get the min and max points
		//CoordBBox boundbox;
		//Coord min, max;
		//boundbox = grid->evalActiveVoxelBoundingBox();
		//min = boundbox.min();
		//max = boundbox.max();
		
		//Force values for zero-centered 1x1x1 cube
		Vec3d min, max;
		min = Vec3d(-0.5,-0.5,-0.5);
		max = Vec3d(0.5,0.5,0.5);		

		//cout << "Min: " << min.x() << " " << min.y() << " " << min.z() << "\n";
	
		int xDim, yDim, zDim;
		float bucketSize = r / sqrt(n);
		xDim = (int)ceil(abs(max.x() - min.x()) / bucketSize); //figure out how many buckets we need of size r/sqrt(n)
		yDim = (int)ceil(abs(max.y() - min.y()) / bucketSize);
		zDim = (int)ceil(abs(max.z() - min.z()) / bucketSize);

		cout << "xDim: " << xDim << " " << "yDim: " << yDim << " " << "zDim: " << zDim << "\n";

		//Init 3d vector of ints and resize each vector to be the right size, then fill with "-1"
		//Height = yDim
		//Width = xDim
		//Depth = zDim (follows Houdini and Unity convention)
		vector<vector<vector<int>>> backGrid;
		backGrid.resize(yDim);
		for(int i = 0; i < yDim; i++){
			backGrid[i].resize(xDim);
			for(int j = 0; j <  xDim; j++){
				backGrid[i][j].resize(zDim);
				for(int k = 0; k < zDim; k++){
					backGrid[i][j][k] = -1;
				}
			}
		}		

		//Initialize a grid sampler so we can find the values of continuous points sampled in the grid!		
		tools::GridSampler<FloatTree, tools::BoxSampler> interpolator(grid->constTree(), grid->transform());
		
		cout << "Begin search for x_0.\n";

		//FIND x0: Randomly search the boundbox until we find a point inside the mesh! (value is negative)
		vector<float> randPoint;
		Vec3d p;
		while(1){
			randPoint = randPointInBBox(min, max);
			p = Vec3d(randPoint[0], randPoint[1], randPoint[2]);
			//if we find a value less than 0, we found a point inside the mesh, x_0! Otherwise keep searching
			if(interpolator.wsSample(p) < 0){
				break;
			}
		}

		cout << "Found x_0!\n";

		//Initialize the sample list (list of actual points we are keeping)
		vector<Vec3d> samples;
		samples.push_back(p); //add x_0 to our samples!

		//Insert x_0 into the background grid
		vector<int> indeces = backgroundGridIndeces(p, min, bucketSize);
		backGrid[indeces[0]][indeces[1]][indeces[2]] = 0; //add index 0 to the backGrid

		//Initialize the "active" list (an array of sample indeces) with this index (0)
		vector<int> activeList; //this list will hold active indeces from the samples vector
		activeList.push_back(0); 

		cout << "Begin sampling.\n";

		//While(active list not empty):
		int randIndex;
		int foundNewSample = 0;
		int x_i;
		Vec3d p_i;
		Vec3d currSample;
		vector<int> currIndex;
		int validSample = 0;
		double xLen, yLen, zLen, xUpdate, yUpdate, zUpdate;
		int lessX, lessY, lessZ, moreX, moreY, moreZ = 0;
		while(!activeList.empty()){

			foundNewSample = 0;

			//Choose random index from the list
			randIndex = 0 + (rand() % static_cast<int>(activeList.size() - 1 - 0 + 1)); //gets rand index between 0 and activeList.size - 1

			x_i = activeList[randIndex]; //index of x_i in sample list
			p_i = samples[x_i]; //actual point we are querying

			//Swap the selected point in activeList with the final element so we can just pop_back if it's a bad point!
			activeList[randIndex] = activeList[activeList.size()-1];
			activeList[activeList.size()-1] = x_i;

			cout << "Begin generating annulus points\n";

			//Generate UP TO K points chosen uniformly from the spherical annulus between radius r and 2r around x_i
			for(int i = 0; i < kParam; i++){

				//Choose a point in the spherical annulus between radius r and 2r around x_i
				//Then we need to move it inside the mesh if it's outside
				currSample = getPointInAnnulus(p_i, r);
				if(interpolator.wsSample(currSample) >= 0){ //if currSample gives a positive value we know it's NOT in the mesh
					cout << "Correcting annulus point\n";
					//cout << "Sample Before Correction: " << currSample.x() << " " << currSample.y() << " "<< currSample.z() << "\n";
					xLen,yLen,zLen;
					xLen = max.x() - min.x();
					yLen = max.y() - min.y();
					zLen = max.z() - min.z();
					if(currSample.x() < min.x()){ lessX = 1; } else{lessX = 0;} //determine which dimensions need to be wrapped around
					if(currSample.y() < min.y()){ lessY = 1; } else{lessY = 0;}
					if(currSample.z() < min.z()){ lessZ = 1; } else{lessZ = 0;}
					if(currSample.x() > max.x()){ moreX = 1; } else{moreX = 0;}
					if(currSample.y() > max.y()){ moreY = 1; } else{moreY = 0;}
					if(currSample.z() > max.z()){ moreZ = 1; } else{moreZ = 0;}
					xUpdate = currSample.x() + (lessX * xLen) - (moreX * xLen);
					yUpdate = currSample.y() + (lessY * yLen) - (moreY * yLen);
					zUpdate = currSample.z() + (lessZ * zLen) - (moreZ * zLen);
					currSample = Vec3d( xUpdate, yUpdate, zUpdate);
					//cout << lessX << " " << lessY << " " << lessZ << " " << moreX << " " << moreY << " " << moreZ << "\n";
					//cout << "Sample After Correction: " << currSample.x() << " " << currSample.y() << " " << currSample.z() << "\n"; 
				}

				if(interpolator.wsSample(currSample) >= 0){
					cout << "Somehow still not inside mesh... :(\n";
					//cout << "Sample After Correction: " << currSample.x() << " " << currSample.y() << " " << currSample.z() << "\n"; 
				}
				
				cout << "Found good annulus point!\n";

				//Check if it is within distance r of existing samples (using background grid to only test neighbors)
				currIndex = backgroundGridIndeces(currSample, min, bucketSize);	//get index of our new possible sample		
				
				cout << "Got currSample's backgrid indeces, now checking neighbor cells\n";
				
				cout << "Backgrid indeces: " << currIndex[0] << " " << currIndex[1] << " " << currIndex[2] << "\n";

				validSample = checkNeighborCells(currIndex, currSample, backGrid, samples, r); //check whether it is within distance r of samples
				
				if(validSample){
					foundNewSample = 1;
					backGrid[currIndex[0]][currIndex[1]][currIndex[2]] = samples.size(); //add new index to backgrid
					activeList.push_back(samples.size()); //add this index to the active list as well!
					samples.push_back(currSample); //add new sample to the sample list (index = size of sample list before adding)
					//cout << "Found a new sample!\n";
					break; //now that we found a sample from the current point we can quit searching for one!
				}	
			}
			//If after K attempts, no such point is found, instead remove x_i from the active list
			if(foundNewSample == 0){
				//cout << "Didn't find a new sample :(\n";
				activeList.pop_back(); //back element is the one to pop since we swapped earlier!
			}
			cout << activeList.size() << "\n";
			
		}

		//Save the points from the samples list 
		cout << "Number of samples: " << samples.size() << "\n";
		saveSamples(samples);
		
	}

	return 0;
}

//Return a random point in the bound box defined by the given min and max points
vector<float> randPointInBBox(Vec3d min, Vec3d max){

	vector<float> point(3,0);

	point[0] = min.x() + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(max.x()-min.x())));
	point[1] = min.y() + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(max.y()-min.y())));
	point[2] = min.z() + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(max.z()-min.z())));

	return point;

}

//Return a vector (i,j,k) that represents the indeces to place "point" into the background grid
//Remember that it indexes y, then x, then z!
vector<int> backgroundGridIndeces(Vec3d point, Vec3d min, float bucketSize){

	vector<int> indeces(3,-1);
	int i, j, k = 0;
	float temp = 0;
	//Quicker way to get grid index!
	i = std::floor((point.y() - min.y())/bucketSize);
	j = std::floor((point.x() - min.x())/bucketSize);
	k = std::floor((point.z() - min.z())/bucketSize);
	/*	
	while(1){
		//Figure out y index first (i)
		temp = min.y() + (bucketSize * i);
		if(temp > point[1]){
			i = i - 1; //we found the bucket! the correct index is one less than the one required to make temp > p.y
			break;
		}
		i = i + 1; //increment i otherwise
	}
	while(1){
		//Figure out x index (j)
		temp = min.x() + (bucketSize * j);
		if(temp > point[0]){
			j = j - 1; //we found the bucket! the correct index is one less than the one required to make temp > p.x
			break;
		}
		j = j + 1; //increment i otherwise
	}
	while(1){
		//Figure out z index finally (k)
		temp = min.z() + (bucketSize * k);
		if(temp > point[2]){
			k = k - 1; //we found the bucket! the correct index is one less than the one required to make temp > p.z
			break;
		}
		k = k + 1; //increment k otherwise
	}
	*/
	indeces[0] = i;
	indeces[1] = j;
	indeces[2] = k;

	return indeces;
}

//Return a point in the spherical annulus around "center" between radius r and 2r
Vec3d getPointInAnnulus(Vec3d center, float r){
	
	//randomly generate spherical coordinates (rho,theta,phi)
	float rho, theta, phi;
	rho = r + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/((2*r) - r))); //between r and 2r
	theta = 0 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2*M_PI - 0))); //between 0 and 2pi
	phi = 0 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2*M_PI - 0))); //between 0 and 2pi

	//Now convert to cartesian coordinates
	//x = rho * sin phi * cos theta
	//y = rho * sin phi * sin theta
	//z = rho * cos phi
	//NOTE: Must also add the coords of the center to translate it in space
	Vec3d point = Vec3d(0,0,0);
	point[0] = center[0] + (rho * sin(phi) * cos(theta));
	point[1] = center[1] + (rho * sin(phi) * sin(theta));
	point[2] = center[2] + rho * cos(phi);
	
	//either add or subtract (max - min) based oin whether we're greater than max or less than min!
	return point;
}

//Check all neighbor cells for existing points and check whether they are sufficiently far apart to exist near each other or not
int checkNeighborCells(vector<int> currIndex, Vec3d currPoint, vector<vector<vector<int>>> backGrid, vector<Vec3d> samples, float r){

	int validSample = 1;
	
	int i, j, k;
	i = currIndex[0];
	j = currIndex[1];
	k = currIndex[2];

	int iLength = backGrid.size();
	int jLength = backGrid[0].size();
	int kLength = backGrid[0][0].size();

	int iIndex, jIndex, kIndex;

	if(backGrid[i][j][k] != -1){
		//if there is already a sample in this cell we instantly know they cannot be far enough apart
		validSample = 0;
		return validSample;
	}

	int point2Ind; //hold index of second neighbor point

	for(int r = i-2; r < i+3; r++){ //iterate through i-2, i-1, i, i+1, and i+2
		for(int s = j-2; s < j+3; s++){ //iterate through j-2, j-1, j, j+1, and j+2
			for(int t = k-2; t < k+3; t++){ //iterate through k-2, k-1, k, k+1, and k+2
				
				//This will set the index to wrap around if we need to wrap around in the positive direction
				iIndex = r % iLength;
				jIndex = s % jLength;
				kIndex = t % kLength;

				//If at an edge, wrap around in negative direction!
				if( r == -1 ){ iIndex = (iLength - 1); }
				if( s == -1 ){ jIndex = (jLength - 1); } 	
				if( t == -1 ){ kIndex = (kLength - 1); }
				if( r == -2 ){ iIndex = (iLength - 2); }
				if( s == -2 ){ jIndex = (jLength - 2); } 	
				if( t == -2 ){ kIndex = (kLength - 2); }		

				cout << "i j k: " << iIndex << " " << jIndex << " " << kIndex << "\n";

				if(backGrid[iIndex][jIndex][kIndex] != -1){
					//if we find a non negative value, we need to check the distance between the two
					point2Ind = backGrid[iIndex][jIndex][kIndex];
					if(!distOk(currPoint, samples[point2Ind], r)){
						//if they are not far enough apart, throw out this sample
						validSample = 0;
						return validSample;
					}
				}
			}	
		}
	}

	//if we make it through to here, then the sample should be good!
	validSample = 1;
	return validSample;

}

int distOk(Vec3d p1, Vec3d p2, float r){

	int farEnough;	
	float dist = sqrt( (p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2]) );
	if(dist <= r){
		farEnough = 0;
	}
	else{
		farEnough = 1; //far enough!
	}

}

void saveSamples(vector<Vec3d> samples){
	
	//Write to a binary file and a text file
	ofstream outfile;
	ofstream outfile2;
	outfile.open("samples.dat", ios::binary | ios::out);
	outfile2.open("samples.txt");
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
	io::File file("points.vdb");
	GridPtrVec grids;
	grids.push_back(grid);
	file.write(grids);
	file.close();	


	return;
}

#pragma once

#include <maya/MPoint.h>
#include <maya/MDagPath.h>
#include <vector>

using std::vector;

class Voxel
{
public:
	Voxel(void);
	~Voxel(void);

	MPoint min;
	MPoint max;
	
	vector<MDagPath> meshes;


};


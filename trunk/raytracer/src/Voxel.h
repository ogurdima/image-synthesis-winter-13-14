#pragma once

#include <maya/MPoint.h>
#include <maya/MDagPath.h>
#include <vector>
#include "Util.h"

using std::vector;
using namespace util;

class Voxel
{
public:
	Voxel(MPoint _min, MPoint _max);
	~Voxel(void);

	MPoint		min;
	MPoint		max;
	vector<int> meshIdx;

	bool intersectsWith(MPoint otherMin, MPoint otherMax);
	bool intersectsWith(MDagPath meshPath, MPoint& intersectionPoint);
};


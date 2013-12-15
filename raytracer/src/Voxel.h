#pragma once

#include <maya/MPoint.h>
#include <maya/MDagPath.h>
#include <maya/MItMeshPolygon.h>
#include <vector>
#include "Plane.h"
#include "Profiler.h"
#include "Util.h"
#include "Mesh.h"

using std::vector;
using namespace util;

class Voxel
{
	MPoint		min;
	MPoint		max;
	MPoint		center;

	vector<Plane> planes;

public:
	Voxel(MPoint _min, MPoint _max);
	~Voxel(void);

	inline MPoint Min()
	{
		return min;
	}

	inline MPoint Max()
	{
		return max;
	}

	bool		intersectsWith(const MPoint& otherMin, const MPoint& otherMax);

	bool		intersectsWith(const MeshDataT& meshPath,const double halfsSides[3], vector<int>& faceIds) const;
	
	bool		findExitDirection(const MPoint& src, const MVector& dir, AxisDirection& farDir);
};


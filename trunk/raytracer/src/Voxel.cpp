#include "Voxel.h"


Voxel::Voxel(MPoint _min, MPoint _max)
{
	min = _min;
	max = _max;
	meshIdx.clear();
}


bool Voxel::intersectsWith(MPoint otherMin, MPoint otherMax)
{
	if (	intervalsOverlap(min.x, max.x, otherMin.x, otherMax.x) && 
			intervalsOverlap(min.y, max.y, otherMin.y, otherMax.y) && 
			intervalsOverlap(min.z, max.z, otherMin.z, otherMax.z))
	{
		return true;
	}
	return false;
}

// TODO: Implement this
bool Voxel::intersectsWith(MDagPath meshPath, MPoint& intersectionPoint)
{
	return false;
}

Voxel::~Voxel(void)
{
};



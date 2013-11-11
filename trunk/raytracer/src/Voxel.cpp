#include "Voxel.h"


Voxel::Voxel(MPoint _min, MPoint _max):
xnext(NULL),
xprev(NULL),
ynext(NULL),
yprev(NULL),
znext(NULL),
zprev(NULL)
{
	min = _min;
	max = _max;
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



#include "Voxel.h"
#include "Plane.h"


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

bool Voxel::intersectionsWithRay( const MPoint& src, const MVector& dirVec, MPoint & nearInt, AxisDirection& nearDir, MPoint farInt, AxisDirection& farDir )
{
	vector<Plane> planes;
	planes.resize(6);

	planes[X_NEG] = Plane(min, MVector(-1,0,0));
	planes[X_POS] = Plane(max, MVector(1,0,0));
	planes[Y_NEG] = Plane(min, MVector(0, -1 ,0));
	planes[Y_POS] = Plane(max, MVector(0,1,0));
	planes[Z_NEG] = Plane(min, MVector(0,0, -1));
	planes[Z_POS] = Plane(max, MVector(0,0,1));

	double times[6];
	AxisDirection dirs[6];
	MPoint ints[6];
	int count = 0;

	for (int i = 0; i < UNKNOWN_DIR && count < 2; ++i)
	{
		dirs[count] = (AxisDirection) i;
		if( planes[i].rayIntersection(src, dirVec, times[count], ints[count]) 
			&& pointInRectangle(dirs[count], ints[count],  min , max))
		{
			 count++;
		}
	}
	if(count != 2)
	{
		return false;
	}
	if(times[0] < times[1])
	{
		nearInt = ints[0];
		nearDir = dirs[0];
		farInt = ints[1];
		farDir = dirs[1];
	}
	else
	{
		nearInt = ints[1];
		nearDir = dirs[1];
		farInt = ints[0];
		farDir = dirs[0];
	}
	return true;
}



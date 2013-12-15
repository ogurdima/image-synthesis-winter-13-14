#include "Voxel.h"

Voxel::Voxel(MPoint _min, MPoint _max)
{
	min = _min;
	max = _max;
	center = (min + max) / 2;
	planes.resize(6);

	planes[X_NEG] = Plane(min, MVector(-1,0,0));
	planes[X_POS] = Plane(max, MVector(1,0,0));
	planes[Y_NEG] = Plane(min, MVector(0, -1 ,0));
	planes[Y_POS] = Plane(max, MVector(0,1,0));
	planes[Z_NEG] = Plane(min, MVector(0,0, -1));
	planes[Z_POS] = Plane(max, MVector(0,0,1));

}


bool Voxel::intersectsWith(const MPoint& otherMin, const MPoint& otherMax)
{
	if (	intervalsOverlap(min.x, max.x, otherMin.x, otherMax.x) && 
			intervalsOverlap(min.y, max.y, otherMin.y, otherMax.y) && 
			intervalsOverlap(min.z, max.z, otherMin.z, otherMax.z))
	{
		return true;
	}
	return false;
}


bool Voxel::intersectsWith(const MeshDataT& mesh, const double halfsSides[3], vector<int>& faceIds) const
{
	int size = mesh.faces.size();
	for (int i = 0; i < size; ++i)
	{
		if(triangleBoxOverlap(center, halfsSides, mesh.faces[i].vertices))
		{
			faceIds.push_back(i);
		}
	}
	
	return !faceIds.empty();
}

Voxel::~Voxel(void)
{
};

bool Voxel::findExitDirection( const MPoint& src, const MVector& dirVec, AxisDirection& farDir )
{
	double times[2];
	AxisDirection dirs[2];
	MPoint intersection;
	int count = 0;

	for (int i = 0; i < UNKNOWN_DIR && count < 2; ++i)
	{
		dirs[count] = (AxisDirection) i;
		if( planes[i].rayIntersection(src, dirVec, times[count], intersection) 
			&& pointInRectangle(dirs[count], intersection,  min , max))
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
		farDir = dirs[1];
	}
	else
	{
		farDir = dirs[0];
	}

	return true;
}



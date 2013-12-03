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


bool Voxel::intersectsWith(const MDagPath& meshPath, const double halfsSides[3], vector<int>& faceIds) const
{
	
	MItMeshPolygon faceIt(meshPath);

	MPointArray ptsArray;
	MIntArray vertexList;
	MPoint pts[3];
	for (;!faceIt.isDone();faceIt.next())
	{
		faceIt.getTriangles(ptsArray,vertexList,MSpace::kWorld);
		if(triangleBoxOverlap(center, halfsSides, ptsArray))
		{
			faceIds.push_back(faceIt.index());
		}
	}
	
	return !faceIds.empty();
}

Voxel::~Voxel(void)
{
};

bool Voxel::intersectionsWithRay( const MPoint& src, const MVector& dirVec, MPoint & nearInt, AxisDirection& nearDir, MPoint farInt, AxisDirection& farDir )
{
	//Profiler::startTimer("SELF::intersectionsWithRay");
	
	double times[2];
	AxisDirection dirs[2];
	MPoint ints[2];
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
		//Profiler::finishTimer("SELF::intersectionsWithRay");
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
	//Profiler::finishTimer("SELF::intersectionsWithRay");
	return true;
}



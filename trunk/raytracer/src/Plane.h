#pragma once

#include <maya/MPoint.h>
#include <maya/MVector.h>
#include "Util.h"

struct Plane {
	
	MPoint planePoint;
	MVector normal;
	
	Plane()
	{
		normal = MVector(1,0,0);
	}

	Plane(const MPoint& _point , const MVector& _normal)
	{
		planePoint = _point;
		normal = _normal;
		normal.normalize();
	}

	bool rayIntersection(const MPoint& rayOrigin, const MVector &rayDirection, double &time, MPoint &intersectionPoint);
};


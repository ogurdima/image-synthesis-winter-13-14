#include "Plane.h"



bool Plane::rayIntersection( const MPoint& rayOrigin, const MVector &rayDirection, double &time, MPoint &intersectionPoint )
{
	double denom = normal *  rayDirection;
	if (abs(denom) > DOUBLE_NUMERICAL_THRESHHOLD) {
		MVector originToPlanePoint =  planePoint - rayOrigin;
		time = (originToPlanePoint * normal) / denom;
		intersectionPoint = rayOrigin + (rayDirection * time );
		return true;
	}
	return false;
}

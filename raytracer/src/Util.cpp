#include "Util.h"

namespace util 
{
	MString pointToString(MPoint p)
	{
		MString outstr;
		outstr += "["; 
		outstr += p.x;
		outstr += ",";
		outstr += p.y;
		outstr += ",";
		outstr += p.z;
		outstr += ",";
		outstr += p.w;
		outstr += "]";
		return outstr;
	}

	MString vectorToString(MVector p)
	{
		MString outstr;
		outstr += "["; 
		outstr += p.x;
		outstr += ",";
		outstr += p.y;
		outstr += ",";
		outstr += p.z;
		outstr += "]";
		return outstr;
	}

	MString colorToString(MColor c)
	{
		MString outstr;
		outstr += "["; 
		outstr += c.r;
		outstr += ",";
		outstr += c.g;
		outstr += ",";
		outstr += c.b;
		outstr += ",";
		outstr += c.a;
		outstr += "]";
		return outstr;
	}

	MMatrix getDagPathTransformationMatrix(MDagPath dagPath, MStatus* statusPtr)
	{
		MStatus status;
		MObject transformNode = dagPath.transform(&status);
		if (status != MS::kSuccess) {
			if (statusPtr != NULL) {
				*statusPtr = status;
			}
			return MMatrix();
		}
		MFnDagNode transform(transformNode, &status);
		if (status != MS::kSuccess) {
			if (statusPtr != NULL) {
				*statusPtr = status;
			}
			return MMatrix();
		}
		MTransformationMatrix matrix(transform.transformationMatrix(&status));
		if (status != MS::kSuccess) {
			if (statusPtr != NULL) {
				*statusPtr = status;
			}
			return MMatrix();
		}
		return matrix.asMatrix();
	}

	pair<MPoint, MPoint> computeWfAxisAlignedBoundingBox(MDagPath meshPath, MStatus* statusPtr)
	{
		MStatus status;
		MFnMesh mesh(meshPath, &status);
		if ( status != MS::kSuccess) {
			if (statusPtr != NULL) {
				*statusPtr = status;
			}
			return pair<MPoint,MPoint>(MPoint(), MPoint());
		}

		MPointArray pts;
		MPoint min = MPoint( DBL_MAX ,DBL_MAX, DBL_MAX);
		MPoint max = MPoint( DBL_MIN ,DBL_MIN, DBL_MIN);
		mesh.getPoints(pts, MSpace::kWorld);
		for(int pi = pts.length() - 1; pi >= 0; --pi)
		{
			MPoint current = pts[pi];
			minimize(&(min.x), current.x);
			minimize(&(min.y), current.y);
			minimize(&(min.z), current.z);

			maximize(&(max.x), current.x);
			maximize(&(max.y), current.y);
			maximize(&(max.z), current.z);
		}
		return pair<MPoint,MPoint>(min, max);
	}

	double minimize(double* oldPtr, double newVal)
	{
		if (oldPtr == NULL) {
			return DBL_MAX;
		}
		if (*oldPtr > newVal) {
			*oldPtr = newVal;
		}
		return *oldPtr;
	}

	double maximize(double* oldPtr, double newVal)
	{
		if (oldPtr == NULL) {
			return DBL_MIN;
		}
		if (*oldPtr < newVal) {
			*oldPtr = newVal;
		}
		return *oldPtr;
	}

	bool intervalsOverlap(double x1, double y1, double x2, double y2)
	{
		double l1, r1, l2, r2;
		if (x1 < y1) {
			l1 = x1;
			r1 = y1;
		}
		else {
			l1 = y1;
			r1 = x1;
		}

		if (x2 < y2) {
			l2 = x2;
			r2 = y2;
		}
		else {
			l2 = y2;
			r2 = x2;
		}

		if (l1 < l2) {
			if (r1 < l2) {
				return false;
			}
			return true;
		}

		if (l2 < l1) {
			if (r2 < l1) {
				return false;
			}
			return true;
		}

	}
}
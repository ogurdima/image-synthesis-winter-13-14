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
		MPoint max = MPoint( -DBL_MAX ,-DBL_MAX, -DBL_MAX);
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
			return -DBL_MAX;
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


		if (r1 < (l2 - DOUBLE_NUMERICAL_THRESHHOLD) || r2 < (l1 - DOUBLE_NUMERICAL_THRESHHOLD)) {
			return false;
		}
		return true;			
	}

	int	flatten3dCubeIndex(int dimSize, int x, int y, int z)
	{
		return x + dimSize*y + dimSize*dimSize*z;
	}

	bool isPointInVolume(const MPoint& point, const MPoint& minVolume, const MPoint& maxVolume)
	{
		for(int i = 0; i < 3; ++i)
		{
			if( point[i] < minVolume[i] || point[i] > maxVolume[i])
				return false;
		}
		return true;
	}

	bool valueInInterval( double value, double intervalMin, double intervalMax )
	{
		return !( value < intervalMin || value > intervalMax);
	}

	bool pointInRectangle( AxisDirection projectionDirection, const MPoint point, const MPoint minPoint, const MPoint maxPoint )
	{
		switch (projectionDirection)
		{
		case X_NEG:
		case X_POS:
			return valueInInterval(point.y, minPoint.y, maxPoint.y) && valueInInterval(point.z, minPoint.z, maxPoint.z);
		case Y_NEG:
		case Y_POS:
			return valueInInterval(point.x, minPoint.x, maxPoint.x) && valueInInterval(point.z, minPoint.z, maxPoint.z);
		case Z_NEG:
		case Z_POS:
			return valueInInterval(point.y, minPoint.y, maxPoint.y) && valueInInterval(point.x, minPoint.x, maxPoint.x);
		default:
			break;
		}
		return false;
	}

}

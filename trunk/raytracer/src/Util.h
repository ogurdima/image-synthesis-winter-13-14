#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <math.h>
#include <maya/MArgList.h>
#include <maya/MObject.h>
#include <maya/MGlobal.h>
#include <maya/MFnMesh.h>
#include <maya/MMatrix.h>
#include <maya/MFnTransform.h>
#include <maya/MPointArray.h>
#include <maya/MDagPath.h>
#include <maya/MFloatMatrix.h>
#include <maya/MSyntax.h>
#include <maya/MFnLight.h>
#include <maya/M3dView.h>
#include <maya/MImage.h>


using std::vector;
using std::string;
using std::ostringstream;
using std::pair;

namespace util
{
	MString							pointToString(MPoint p);
	MString							vectorToString(MVector p);
	MString							colorToString(MColor c);
	MMatrix							getDagPathTransformationMatrix(MDagPath dagPath, MStatus* statusPtr);
	pair<MPoint, MPoint>			computeWfAxisAlignedBoundingBox(MDagPath meshPath, MStatus* statusPtr = NULL);
	inline double					minimize(double* oldPtr, double newVal);
	inline double					maximize(double* oldPtr, double newVal);

	bool							intervalsOverlap(double x1, double y1, double x2, double y2);

};


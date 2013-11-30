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
#include <maya/MFnLambertShader.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include "Definitions.h"
#include "Profiler.h"

using std::vector;
using std::string;
using std::ostringstream;
using std::pair;

const double					DOUBLE_NUMERICAL_THRESHHOLD = 0.0000001;

namespace util
{
	

	MString							pointToString(MPoint p);
	MString							vectorToString(MVector p);
	MString							colorToString(MColor c);
	//MMatrix							getDagPathTransformationMatrix(MDagPath dagPath, MStatus* statusPtr);
	pair<MPoint, MPoint>			computeWfAxisAlignedBoundingBox(MDagPath meshPath, MStatus* statusPtr = NULL);
	inline double					minimize(double* oldPtr, double newVal);
	inline double					maximize(double* oldPtr, double newVal);

	bool							valueInInterval(double value, double intervalMin, double intervalMax);
	bool							intervalsOverlap(double x1, double y1, double x2, double y2);
	int								flatten3dCubeIndex(int dimSize, int x, int y, int z);
	bool							pointInRectangle(AxisDirection projectionDirection, const MPoint point, const MPoint minPoint, const MPoint maxPoint );
	bool							isPointInVolume(const MPoint& point, const MPoint& minVolume, const MPoint& maxVolume);
	bool							triangleBoxOverlap( const MPoint center , const double boxhalfsize[3], const MPointArray triangleVertices);
	bool							rayIntersectsTriangle(const MPoint raySrc,const MVector rayDirection, const MPoint triangleVertices[3], double& time, MPoint& intersection);
	MVector							reflectedRay(MVector ligthDir, MVector normal);
	void							caclulateBaricentricCoordinates( MPoint triangleVertices[3], MPoint point, double baricentricCoords[3] );

	bool							getLambertShaderTexture(MFnLambertShader& lambert, MImage& img);

	MColor							sumColors(const MColor& c1 , const MColor& c2);
	MColor							textureNearesNeighborAtPoint(const MImage* texture, double u, double v, bool repeat = true);

};


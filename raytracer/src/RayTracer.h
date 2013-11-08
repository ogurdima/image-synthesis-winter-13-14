#pragma once

#include <maya/MArgList.h>
#include <maya/MObject.h>
#include <maya/MGlobal.h>
#include <maya/MPxCommand.h>
#include <maya/MItDependencyNodes.h>
#include <maya/MItGeometry.h>
#include <maya/MFnMesh.h>
#include <maya/MItMeshVertex.h>
#include <maya/MMatrix.h>
#include <maya/MFnTransform.h>
#include <maya/MFnCamera.h>
#include <maya/MPointArray.h>
#include <maya/MItDag.h>
#include <maya/MDagPath.h>
#include <maya/MFloatMatrix.h>
#include <maya/MSyntax.h>
#include <maya/MFnLight.h>
#include <maya/M3dView.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatVector.h>
#include <maya/MImage.h>



class RayTracer : public MPxCommand 
{

	int width;
	int height;

	MVector eyePosition;

	MFloatVector ** rayDirections;

public:

	RayTracer() {
		width = 640;
		height = 480;
	};
	virtual MStatus doIt(const MArgList& argList);
	static void* creator();

	void printMeshPoints();
	void printCamerasInfo();
	void printObjectTypesInScene2();
	void printObjectTypesInScene();
	MPoint getObjectSpaceCentroid(MObject obj);
	MString pointToStr(MPoint p);
	

	void getCameraInfo();
	void goOverRays();

};
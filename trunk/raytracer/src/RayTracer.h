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
#include <maya/MBoundingBox.h>
#include <maya/MSelectionList.h>
#include <vector>
#include <string>
#include <sstream>

#include "Util.h"


using std::vector;
using std::string;
using std::ostringstream;
using std::pair;

using namespace util;

class RayTracer : public MPxCommand
{

	MPoint minScene;
	MPoint maxScene;

	MVector eyePosition;
	MFloatVector ** rayDirections;


	int imgWidth;
	int imgHeight;


	struct CameraDataT
	{
		MPoint		eye;
		MVector		viewDir;
		MVector		upDir;
		double		focalLengthCm;
		double		filmWidthCm;

		MString		toString() {
			ostringstream os;
			os << "Camera: Eye:" << pointToString(eye) << ",View:" << vectorToString(viewDir) << ",Up:" 
				<< vectorToString(upDir) << ",FocalDist:" << focalLengthCm << ",FilmWidth:" << filmWidthCm;
			return MString(os.str().c_str());
		}

	} activeCameraData;

	struct ImagePlaneDataT
	{
		MVector		x;
		MVector		y;
		MPoint		lt; //left top
		MPoint		rt; //right top
		MPoint		lb; //left bottom
		MPoint		rb; //right bottom
		double		dp; //delta p - pixel size
	} imagePlane;

	struct LightDataT
	{
		enum { UNDEF, AMBIENT, DIRECTIONAL, POINT} type;
		MPoint		position;
		MVector		direction;
		MColor		color;

		LightDataT() : type(UNDEF) {}
		MString		toString() {
			ostringstream os;
			string typeStr = (type==AMBIENT) ? "Ambient" : (type==DIRECTIONAL) ? "Directional" : (type==POINT) ? "Point" : "Undef";
			os << "Type:" << typeStr;
			if (type != UNDEF) {
				os << ",Col:" << colorToString(color);
			}
			if (type == DIRECTIONAL) {
				os << ",Dir:" << vectorToString(direction);
			}
			if (type == POINT) {
				os << ",Pos:" << pointToString(position);
			}
			return MString(os.str().c_str());
		}
		
	};

	vector<LightDataT> lightingData;

	struct MeshDataT
	{
		MDagPath	dagPath;
		MPoint		max;		// WS axis aligned bounding box min
		MPoint		min;		// WS axis aligned bounding box max

	};

	vector<MeshDataT> meshesData;

	struct VoxelezationParamT
	{
		double		dx;
		double		dy;
		double		dz;
		int			voxelsPerDimension;

		VoxelezationParamT() : voxelsPerDimension(1), dx(1), dy(1), dz(1) {}

	} voxelParams;


public:

	RayTracer() {
		imgWidth  = 640;
		imgHeight = 480;
		minScene = MPoint( DBL_MAX ,DBL_MAX,DBL_MAX);
		maxScene = MPoint(DBL_MIN, DBL_MIN, DBL_MIN);
		voxelParams.voxelsPerDimension = 5;

		meshesData.clear();
		lightingData.clear();

	};
	virtual MStatus doIt(const MArgList& argList);
	static void* creator();

	void printMeshPoints();
	void printCamerasInfo();
	void printObjectTypesInScene2();
	void printObjectTypesInScene();
	MPoint getObjectSpaceCentroid(MObject obj);


	void getCameraInfo();
	void goOverRays();
	void calculateSceneBoundingBox();
	



	void storeActiveCameraData();
	void computeAndStoreImagePlaneData();
	void storeLightingData();
	void storeAmbientLight(MDagPath lightDagPath);
	void storeDirectionalLight(MDagPath lightDagPath);
	void storePointLight(MDagPath lightDagPath);
	void computeAndStoreMeshData();
	void computeAndStoreSceneBoundingBox();
	void voxelizeScene();
	void computeAndStoreVoxelParams();


};
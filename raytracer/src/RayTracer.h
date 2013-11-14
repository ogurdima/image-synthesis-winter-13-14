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
#include <maya/MTimer.h>
#include <vector>
#include <string>
#include <sstream>
#include "Plane.h"
#include "Definitions.h"
#include <map>

#include "Util.h"
#include "Voxel.h"


using std::vector;
using std::string;
using std::ostringstream;
using std::pair;
using std::map;

using namespace util;

class RayTracer : public MPxCommand
{



	MPoint minScene;
	MPoint maxScene;
	vector<Plane> sceneBBPlanes;


	int imgWidth;
	int imgHeight;

	bool cameraInSceneBB;
	int initCameraVoxelX;
	int initCameraVoxelY;
	int initCameraVoxelZ;


	

	

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
		//double		dx;
		//double		dy;
		//double		dz;

		double dimensionDeltaHalfs[3];
		double dimensionDeltas[3];

		int			voxelsPerDimension;

		VoxelezationParamT() : voxelsPerDimension(1)
		{
			dimensionDeltaHalfs[0] = dimensionDeltaHalfs[1] = dimensionDeltaHalfs[2] = 0.5;
			dimensionDeltas[0] = dimensionDeltas[1] = dimensionDeltas[2] = 1;
		}

	} voxelParams;

	struct VoxelDataT
	{
		Voxel* v;
		vector<int> containedMeshIndexes;
		map<int, vector<int>> meshIdToFaceIds;
		
	};

	vector<VoxelDataT> voxelsData;


public:

	RayTracer();
	~RayTracer();
	virtual MStatus doIt(const MArgList& argList);
	static void* creator();

	void triangulateMesh(const MFnMesh& mesh);
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
	void computeAndStoreRawVoxelsData();
	void computeVoxelNeighborhoodData();
	void computeVoxelMeshBboxIntersections();
	void bresenhaim();

	void shootRay(int dimension, int x, int y, int z, MPoint raySource, MVector rayDirection, unsigned char* pixels, int h, int w );

	bool findStartingVoxelIndeces(const MVector& rayDirection, int& bx, int& by, int& bz);

	bool findIndecesByDimension(const MPoint& point, AxisDirection direction,  int& x, int& y, int& z );

	void foo( AxisDirection direction, MPoint closestIntersection, int * bx, int * by, int * bz );
	void initIndeces( AxisDirection direction, int& x, int& y, int& z );
	void orthonormalDirections( AxisDirection direction, AxisDirection& uDirection, AxisDirection& vDirection );
	bool pointInVoxelByDirection( const MPoint& closestIntersection,VoxelDataT voxel, AxisDirection uDirection );
	void incrementIndeces( AxisDirection uDirection, int& x, int& y, int& z );
};
#pragma once

#include <maya/MArgList.h>
#include <maya/MGlobal.h>
#include <maya/MPxCommand.h>
#include <maya/MFnMesh.h>
#include <maya/MMatrix.h>
#include <maya/MFnCamera.h>
#include <maya/MPointArray.h>
#include <maya/MItDag.h>
#include <maya/MDagPath.h>
#include <maya/MSyntax.h>
#include <maya/MFnLight.h>
#include <maya/M3dView.h>
#include <maya/MImage.h>
#include <maya/MSelectionList.h>
#include <maya/MTimer.h>
#include <maya/MArgParser.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MFnLambertShader.h>
#include <maya/MFnReflectShader.h>
#include <maya/MFnBlinnShader.h>
#include <maya/MFnPhongShader.h>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include "Plane.h"
#include "Definitions.h"
#include "Util.h"
#include "Voxel.h"
#include "Mesh.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


using std::vector;
using std::string;
using std::ostringstream;
using std::pair;
using std::map;

using namespace util;


#define		widthFlag				"-w"
#define		heightFlag				"-h"
#define		voxelsFlag				"-n"
#define		supersamplingFlag		"-s"
#define		superSamplingTypeFlag	"-ss"
#define		rayDepthFlag			"-rd"

#define		RAND_PRECISION			1000
#define		RAND					((double)( rand() % RAND_PRECISION)) / ((double) (RAND_PRECISION - 1))

#define		BACKGROUND_COLOR		MColor(0, 0, 0, 1)

class RayTracer : public MPxCommand
{
	MPoint minScene;
	MPoint maxScene;
	vector<Plane> sceneBBPlanes;

	//int imgWidth;
	//int imgHeight;
	bool cameraInSceneBB;
	int initCameraVoxelX;
	int initCameraVoxelY;
	int initCameraVoxelZ;
	//int supersamplingCoeff;

	static char* outputFilePath;
	static char* statisticsFilePath;

	static double prepTime;
	static double totalTime;
	static double timePerPixel;
	static double timePerPixelStandardDeviation;
	static long intersectionTestCount;
	static long intersectionFoundCount; 
	static long voxelsTraversed;
	static long totalRayCount;
	static long totalPolyCount;
	 
	struct CameraDataT
	{
		MPoint		eye;
		MVector		viewDir;
		MVector		upDir;
		double		focalLengthCm;
		double		filmWidthCm;
		bool		isPerspective;

		MString		toString() {
			ostringstream os;
			os << "Camera: Eye:" << pointToString(eye) << ",View:" << vectorToString(viewDir) << ",Up:" 
				<< vectorToString(upDir) << ",FocalDist:" << focalLengthCm << ",FilmWidth:" << filmWidthCm;
			return MString(os.str().c_str());
		}

	};

	struct ImagePlaneDataT
	{
		enum { UNIFORM, JITTERED, RANDOM, ADAPTIVE, UNDEFINED} ssType;

		int			imgWidth;
		int			imgHeight;


		MVector		x;
		MVector		y;

		MPoint		lt; //left top
		MPoint		rt; //right top
		MPoint		lb; //left bottom
		MPoint		rb; //right bottom
		double		dp; //delta p - pixel size

		MVector		dx;
		MVector		dy;


		// Super sampling params
		double		ssDp;
		int			supersamplingCoeff;
		MVector		ssdx;
		MVector		ssdy;
		

		ImagePlaneDataT()
		{
			srand (time(NULL));
		}

		void	getPointsOnIP(const int w, const int h, vector<MPoint>& out ) const
		{
			out.clear();
			switch (ssType)
			{
			case RayTracer::ImagePlaneDataT::UNIFORM: {
					MPoint pixelLB = lb + h * dy + w * dx;
					for(int ww = 1; ww <= supersamplingCoeff; ++ ww) {
						for(int hh = 1; hh <= supersamplingCoeff; ++hh) {
							out.push_back(pixelLB + hh * ssdx + ww * ssdy);
						}
					}
				}
				break;
			case RayTracer::ImagePlaneDataT::JITTERED: {
					MPoint pixelLB = lb + h * dy + w * dx;
					for(int ww = 1; ww <= supersamplingCoeff; ++ ww) {
						for(int hh = 1; hh <= supersamplingCoeff; ++hh) {
							out.push_back(pixelLB + ((double) hh + RAND) * ssdx + ((double) ww + RAND) * ssdy);
						}
					}
				}

				break;
			case RayTracer::ImagePlaneDataT::RANDOM: {
					MPoint pixelLB = lb + h * dy + w * dx;
					for(int r = 0; r < supersamplingCoeff; ++r) {
						out.push_back(pixelLB + (RAND * dx) + (RAND * dy));
					}
				}
				break;
			case RayTracer::ImagePlaneDataT::ADAPTIVE:
				break;
			default:
				break;
			}
		}

	};

	struct LightDataT
	{
		enum { UNDEF, AMBIENT, DIRECTIONAL, POINT} type;
		MPoint		position;
		MVector		direction;
		MColor		color;
		float		intencity;

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
		
		MVector directionToPoint(const MPoint & point)
		{
			if(type == DIRECTIONAL)
				return direction;
			if(type == POINT)
				return (point - position).normal();
			return MVector(0,0,-1);
		}

		double distanceToPoint(const MPoint & point)
		{
			if(type == POINT)
				return (point - position).length();
			return DBL_MAX;
		}

	};

	

	

	struct SceneParamT
	{
		//double		dx;
		//double		dy;
		//double		dz;

		int rayDepth;

		double dimensionDeltaHalfs[3];
		double dimensionDeltas[3];

		int			voxelsPerDimension;
		int			voxelsPerDimensionSqr;

		SceneParamT() : voxelsPerDimension(1), voxelsPerDimensionSqr(1)
		{
			dimensionDeltaHalfs[0] = dimensionDeltaHalfs[1] = dimensionDeltaHalfs[2] = 0.5;
			dimensionDeltas[0] = dimensionDeltas[1] = dimensionDeltas[2] = 1;
		}

		inline int	flatten3dCubeIndex( int x, int y, int z)
		{
			return x + voxelsPerDimension*y + voxelsPerDimensionSqr*z;
		}

		inline void incrementIndeces( AxisDirection uDirection, int& x, int& y, int& z, int& cur3dIndex )
		{
			switch (uDirection)
			{
			case X_POS:
				++x;
				++cur3dIndex;
				break;
			case X_NEG:
				--x;
				--cur3dIndex;
				break;
			case Y_POS:
				++y;
				cur3dIndex += voxelsPerDimension;
				break;
			case Y_NEG:
				--y;
				cur3dIndex -= voxelsPerDimension;
				break;
			case Z_POS:
				++z;
				cur3dIndex += voxelsPerDimensionSqr;
				break;
			case Z_NEG:
				--z;
				cur3dIndex -= voxelsPerDimensionSqr;
				break;
			default:
				break;
			}
			voxelsTraversed++;
		}

	} ;

	struct VoxelDataT
	{
		Voxel* v;
		//vector<int> containedMeshIndexes;
		map<int, vector<int>> meshIdToFaceIds;
		
	};

	CameraDataT activeCameraData;
	ImagePlaneDataT imagePlane;
	SceneParamT sceneParams;
	vector<MeshDataT> meshesData;
	vector<VoxelDataT> voxelsData;
	vector<LightDataT> lightingData;
public:

#pragma region INTERACTION
	RayTracer();
	~RayTracer();
	virtual MStatus doIt(const MArgList& argList);
	static void* creator();
	static MSyntax newSyntax();
	bool parseArgs( const MArgList& args);
	void openImageInMaya();
	void printStatisticsReport();
#pragma endregion

#pragma region MESH
	void triangulateMesh(const MFnMesh& mesh);
	void computeAndStoreMeshData();
	void computeVoxelMeshIntersections();
	void storeMeshMaterial(MeshDataT& m, const MDagPath& path);
#pragma endregion 

#pragma region LIGHTS
	void storeLightingData();
	void storeAmbientLight(MDagPath lightDagPath);
	void storeDirectionalLight(MDagPath lightDagPath);
	void storePointLight(MDagPath lightDagPath);
#pragma endregion 

#pragma region CAMERA_AND_IMG_PLANE
	void storeActiveCameraData();
	void storeCameraData( MFnCamera &camera );
	void computeAndStoreImagePlaneData();
#pragma endregion 

#pragma region SCENE
	void computeAndStoreSceneBoundingBox();
	void voxelizeScene();
	void computeAndStoreVoxelParams();
	void computeAndStoreRawVoxelsData();
#pragma endregion 

#pragma region ALGO
	void bresenhaim();
	MColor shootRay(const MPoint& raySrc, const MVector& rayDir, int depth);
	bool closestIntersection(const MPoint& raySource,const MVector& rayDirection,int& x,int& y,int& z , int& meshIndex, int& innerFaceId, MPoint& intersection, double depth = DBL_MAX );
	bool closestIntersectionInVoxel(const MPoint& raySource, const MVector& rayDirection, VoxelDataT &voxelData, int &meshIndex, int &innerFaceId, MPoint &intersection);

	bool findStartingVoxelIndeces(const MPoint& raySrc, const MVector& rayDirection, int& bx, int& by, int& bz);

	bool findIndecesByDimension(const MPoint& point, AxisDirection direction,  int& x, int& y, int& z );

	void initIndeces( AxisDirection direction, int& x, int& y, int& z );

	void orthonormalDirections( AxisDirection direction, AxisDirection& uDirection, AxisDirection& vDirection );

	bool pointInVoxelByDirection( const MPoint& closestIntersection,VoxelDataT& voxel, AxisDirection uDirection );

	MColor calculateSpecularAndDiffuse(const MVector& viewDirection,const MVector& lightDirection,const  MVector& normalAtPoint,const MColor& lightColor,const MeshDataT& mesh,const double* bc);
#pragma endregion 

	

};
#include "RayTracer.h"

#include "Profiler.h"
#include "Material.h"

#ifdef _DEBUG
#define DEBUG_REPORT 0
#endif

#ifdef PRINT_FOR_DEBUG
#define PRINT_IN_MAYA(arg) MGlobal::displayInfo((arg));
#endif



#pragma region STATIC VARS
const MString CAMERA_NAME = "cameraShape1";

double	RayTracer::prepTime = 0;
double	RayTracer::totalTime = 0;
double	RayTracer::timePerPixel = 0;
double	RayTracer::timePerPixelStandardDeviation = 0;
long	RayTracer::intersectionTestCount = 0;
long	RayTracer::intersectionFoundCount = 0;
long	RayTracer::voxelsTraversed = 0;
long	RayTracer::totalRayCount = 0;
long	RayTracer::totalPolyCount = 0;
long	RayTracer::totalDepths = 0;
long	RayTracer::totalSamples = 0;
double	RayTracer::samplesPerPixel = 0;
double	RayTracer::samplesPerPixelStdDeviation = 0;




char*	RayTracer::outputFilePath = "C://temp//scene.iff";
char*	RayTracer::statisticsFilePath = "C://temp//stat.txt";

#pragma endregion


MStatus RayTracer::doIt(const MArgList& argList)
{
	cout << "Running raytracer plugin..." << endl;
	MGlobal::displayInfo("Running raytracer plugin...");

	Profiler::clear();
	Profiler::startTimer("doIt::totalTime");
	Profiler::startTimer("doIt::prepTime");

	parseArgs(argList);
	storeActiveCameraData();
	computeAndStoreImagePlaneData();
	storeLightingData();
	computeAndStoreMeshData();
	computeAndStoreSceneBoundingBox();
	voxelizeScene();
	prepTime = Profiler::finishTimer("doIt::prepTime");

	bresenhaim();

	totalTime = Profiler::finishTimer("doIt::totalTime");

	printStatisticsReport();

	openImageInMaya();

	MGlobal::displayInfo("Raytracer plugin run finished!");
	cout << "Raytracer plugin run finished!" << endl;
	return MS::kSuccess;
}

MSyntax RayTracer::newSyntax()
{
	MSyntax syntax;

	syntax.addFlag(widthFlag, "-widthFlag", MSyntax::kLong);
	syntax.addFlag(heightFlag, "-heightFlag", MSyntax::kLong);
	syntax.addFlag(voxelsFlag, "-voxelsFlag", MSyntax::kLong);
	syntax.addFlag(rayDepthFlag, "-rayDepthFlag", MSyntax::kLong);

	syntax.addFlag(supersamplingFlag, "-supersamplingFlag", MSyntax::kLong);
	syntax.addFlag(samplingRateFlag, "-samplingRateFlag", MSyntax::kLong);

	syntax.addFlag(superSamplingTypeFlag, "-superSamplingTypeFlag", MSyntax::kString);
	syntax.addFlag(toleranceFlag, "-toleranceFlag", MSyntax::kDouble);
	syntax.addFlag("-mi", maxSamplingRateFlag, MSyntax::kLong);
	syntax.addFlag("-ma", minSamplingRateFlag, MSyntax::kLong);

	return syntax;
}

bool RayTracer::parseArgs( const MArgList& args)
{
	MArgParser		argData(syntax(), args);
	MStatus			s;

	if ( argData.isFlagSet(widthFlag) ) {
		uint arg;
		s = argData.getFlagArgument(widthFlag, 0, arg);	
		if (s == MStatus::kSuccess) {
			imagePlane.imgWidth = (arg < 1) ? 1 : arg;
		}
	}

	if ( argData.isFlagSet(heightFlag) ) {
		uint arg;
		s = argData.getFlagArgument(heightFlag, 0, arg);	
		if (s == MStatus::kSuccess) {
			imagePlane.imgHeight = (arg < 1) ? 1 : arg;
		}
	}

	if ( argData.isFlagSet(voxelsFlag) ) {
		uint arg;
		s = argData.getFlagArgument(voxelsFlag, 0, arg);	
		if (s == MStatus::kSuccess) {
			sceneParams.voxelsPerDimension = (arg < 1) ? 1 : arg;
			sceneParams.voxelsPerDimensionSqr = sceneParams.voxelsPerDimension * sceneParams.voxelsPerDimension;
		}
	}

	if ( argData.isFlagSet(supersamplingFlag) ) {
		uint arg;
		s = argData.getFlagArgument(supersamplingFlag, 0, arg);	
		if (s == MStatus::kSuccess) {
			imagePlane.supersamplingCoeff = (arg < 1) ? 1 : arg;
		}
	}

	if ( argData.isFlagSet(samplingRateFlag) ) {
		uint arg;
		s = argData.getFlagArgument(samplingRateFlag, 0, arg);	
		if (s == MStatus::kSuccess) {
			imagePlane.supersamplingCoeff = (arg < 1) ? 1 : arg;
		}
	}

	if ( argData.isFlagSet(superSamplingTypeFlag) ) {
		MString arg;
		s = argData.getFlagArgument(superSamplingTypeFlag, 0, arg);	
		if (s == MStatus::kSuccess) {
			if(arg == "uniform")
				imagePlane.ssType = RayTracer::ImagePlaneDataT::UNIFORM;
			else if(arg == "jittered")
				imagePlane.ssType = RayTracer::ImagePlaneDataT::JITTERED;
			else if(arg == "random")
				imagePlane.ssType = RayTracer::ImagePlaneDataT::RANDOM;
			else if(arg == "adaptive")
				imagePlane.ssType = RayTracer::ImagePlaneDataT::ADAPTIVE;
		}
	}

	if(argData.isFlagSet(rayDepthFlag)) {
		uint depth;
		s = argData.getFlagArgument(rayDepthFlag, 0, depth);
		if (s == MStatus::kSuccess && depth >= 1) {
			sceneParams.rayDepth = depth;
		}
	}

	if ( argData.isFlagSet(toleranceFlag) ) {
		double arg;
		s = argData.getFlagArgument(toleranceFlag, 0, arg);	
		if (s == MStatus::kSuccess) {
			imagePlane.ssAdaptiveTolerance = arg;
		}
	}

	if ( argData.isFlagSet(maxSamplingRateFlag) ) {
		uint arg;
		s = argData.getFlagArgument(maxSamplingRateFlag, 0, arg);	
		if (s == MStatus::kSuccess) {
			imagePlane.ssAdaptiveMaxSamples = (arg < 1) ? 1 : arg;
			if (imagePlane.ssAdaptiveMaxSamples > 128) {
				imagePlane.ssAdaptiveMaxSamples = 128;
			}
		}
	}

	if ( argData.isFlagSet(minSamplingRateFlag) ) {
		uint arg;
		s = argData.getFlagArgument(minSamplingRateFlag, 0, arg);	
		if (s == MStatus::kSuccess) {
			imagePlane.ssAdaptiveMinSamples = (arg < 1) ? 1 : arg;
		}
	}

	return true;
}

RayTracer::RayTracer()
{
	imagePlane.imgWidth = 1920;
	imagePlane.imgHeight = 1080;
	imagePlane.supersamplingCoeff = 1;
	imagePlane.ssType = RayTracer::ImagePlaneDataT::UNIFORM;

	imagePlane.ssAdaptiveTolerance = 0.05;
	imagePlane.ssAdaptiveMinSamples = 1;
	imagePlane.ssAdaptiveMaxSamples = 128;
	imagePlane.ssAdaptiveErrorProbability = 0.1; 

	sceneParams.voxelsPerDimension = 1;
	sceneParams.voxelsPerDimensionSqr = 1;
	sceneParams.rayDepth = 1;

	prepTime = 0;
	totalTime = 0;
	timePerPixel = 0;
	timePerPixelStandardDeviation = 0;
	intersectionTestCount = 0;
	intersectionFoundCount = 0;
	voxelsTraversed = 0;
	totalRayCount = 0;
	totalPolyCount = 0;
	totalDepths = 0;
	totalSamples = 0;
	samplesPerPixel = 0;
	samplesPerPixelStdDeviation = 0;

	minScene = MPoint( DBL_MAX ,DBL_MAX,DBL_MAX);
	maxScene = MPoint(-DBL_MAX, -DBL_MAX, -DBL_MAX);

	meshesData.clear();
	lightingData.clear();

	for (int i = 0; i < voxelsData.size(); i++) {
		if (voxelsData[i].v != NULL) {
			delete voxelsData[i].v;
		}
	}
	voxelsData.clear();
};

RayTracer::~RayTracer()
{
	for (int i = 0; i < voxelsData.size(); i++) {
		if (voxelsData[i].v != NULL) {
			delete voxelsData[i].v;
		}
	}
	voxelsData.clear();
}

void* RayTracer::creator()
{
	return new RayTracer;
}

inline void RayTracer::triangulateMesh(const MFnMesh& mesh)
{
	MString cmd("polyTriangulate -ch 0 ");
	cmd += mesh.name();
	MGlobal::executeCommand( cmd );
}

void RayTracer::openImageInMaya()
{
	MString cmd(" file -import -type \"image\" -rpr \"scene\" \"");
	cmd += outputFilePath;
	cmd += "\" ";
	MGlobal::executeCommand( cmd );
}

void RayTracer::printStatisticsReport()
{
	ostringstream os;

	os << "prepTime " << prepTime << endl;
	os << "renderTime " << (totalTime - prepTime) << endl;
	os << "totalTime " << totalTime << endl;
	os << "timePerPixel " << timePerPixel << endl;
	os << "timePerPixelDeviation " << timePerPixelStandardDeviation << endl;
	os << "polygons " << totalPolyCount << endl;
	os << "polygonsPerRay " << ((double)intersectionTestCount / (double)totalRayCount) << endl;
	os << "voxelsPerRay " << (voxelsTraversed / (double)totalRayCount) << endl;
	os << "intersectionTests " << intersectionTestCount << endl;
	os << "Intersections " << ((double)intersectionFoundCount/intersectionTestCount) * 100 << "%" << endl;

	os << "averageSamplingRate " << samplesPerPixel << endl;
	os << "samplingRateDeviation " << samplesPerPixelStdDeviation << endl;
	os << "averageLength "  << totalDepths / (double)totalSamples << endl;

	std::ofstream outfile;
	outfile.open(statisticsFilePath);
	outfile << os.str().c_str();

#ifdef DEBUG_REPORT
	std::cout << os.str().c_str();
#endif

	outfile.close();
}

void RayTracer::storeCameraData( MFnCamera &camera )
{
	activeCameraData.upDir			= camera.upDirection(MSpace::kWorld).normal();
	activeCameraData.viewDir		= camera.viewDirection(MSpace::kWorld).normal();
	activeCameraData.eye			= camera.eyePoint(MSpace::kWorld);
	

	if (camera.isOrtho() ) {
		activeCameraData.isPerspective	= false;
		activeCameraData.filmWidthCm	= camera.orthoWidth();
		activeCameraData.focalLengthCm	= 0.f;
	}
	else {
		activeCameraData.isPerspective	= true;
		activeCameraData.filmWidthCm	= 2.54 * camera.horizontalFilmAperture();
		activeCameraData.focalLengthCm	= (camera.focalLength() / 10);
	}
}

void RayTracer::storeActiveCameraData()
{

	MItDag dagIterator(MItDag::kDepthFirst, MFn::kCamera);
	for(; !dagIterator.isDone(); dagIterator.next())
	{
		MDagPath dagPath;
		dagIterator.getPath(dagPath);
		MFnCamera cur(dagPath);
		if(cur.name() == CAMERA_NAME)
		{
			storeCameraData(cur);
			return;
		}
	}

	MDagPath cameraPath;
	M3dView::active3dView().getCamera( cameraPath );
	MFnCamera camera(cameraPath);
	storeCameraData(camera);

#ifdef PRINT_FOR_DEBUG
	PRINT_IN_MAYA(activeCameraData.toString());
#endif
}

void RayTracer::computeAndStoreImagePlaneData()
{

	imagePlane.x = (activeCameraData.viewDir ^ activeCameraData.upDir).normal();
	imagePlane.y = (imagePlane.x ^ activeCameraData.viewDir).normal();

	MPoint centerPoint = activeCameraData.eye + ( activeCameraData.viewDir * activeCameraData.focalLengthCm );


	double imgAspect = (double)imagePlane.imgWidth / (double)imagePlane.imgHeight;
	//double pixelWidth = (activeCameraData.filmWidthCm)/(double)imgWidth;
	//double pixelHeight = pixelWidth / imgAspect;
	imagePlane.dp = ( activeCameraData.filmWidthCm)/(double)imagePlane.imgWidth;


	imagePlane.ssDp = imagePlane.dp / (double) (imagePlane.supersamplingCoeff + 1);
	imagePlane.ssdx = imagePlane.x * imagePlane.ssDp;
	imagePlane.ssdy = imagePlane.y * imagePlane.ssDp;

	imagePlane.dx = imagePlane.x * imagePlane.dp;
	imagePlane.dy = imagePlane.y * imagePlane.dp;

	int halfWidth = imagePlane.imgWidth / 2;
	int halfHeight = imagePlane.imgHeight / 2;

	imagePlane.lb = centerPoint - (imagePlane.dp * halfWidth * imagePlane.x) - (imagePlane.dp * halfHeight * imagePlane.y);
	imagePlane.lt = centerPoint - (imagePlane.dp * halfWidth * imagePlane.x) + (imagePlane.dp * halfHeight * imagePlane.y);
	imagePlane.rb = centerPoint + (imagePlane.dp * halfWidth * imagePlane.x) - (imagePlane.dp * halfHeight * imagePlane.y);
	imagePlane.rt = centerPoint + (imagePlane.dp * halfWidth * imagePlane.x) + (imagePlane.dp * halfHeight * imagePlane.y);

}

#pragma region storeLighting

void RayTracer::storeLightingData()
{
	MStatus status;
	MItDag dagIterator(MItDag::kDepthFirst, MFn::kLight , &status);
	CHECK_MSTATUS(status);
	for(; !dagIterator.isDone(); dagIterator.next())
	{
		MDagPath dagPath;
		status = dagIterator.getPath(dagPath);
		CHECK_MSTATUS(status);

		if(!dagPath.hasFn(MFn::kLight)) {
			continue;
		}

		if(dagPath.hasFn(MFn::kAmbientLight))
		{
			//PRINT_IN_MAYA("kAmbientLight");
			storeAmbientLight(dagPath);
		}
		else if(dagPath.hasFn(MFn::kDirectionalLight))
		{
			//PRINT_IN_MAYA("kDirectionalLight");
			storeDirectionalLight(dagPath);
		}
		else if(dagPath.hasFn(MFn::kPointLight))
		{
			//PRINT_IN_MAYA("kPointLight");
			storePointLight(dagPath);
		}
		else
		{
#ifdef PRINT_FOR_DEBUG
			PRINT_IN_MAYA("Unsupported light");
#endif
		}
	}
#ifdef PRINT_FOR_DEBUG
	for (int i = 0; i < lightingData.size(); i++)
	{
		PRINT_IN_MAYA(lightingData[i].toString());
	}
#endif
}

void RayTracer::storeAmbientLight(MDagPath lightDagPath)
{
	MStatus status;
	MFnLight l(lightDagPath);
	MColor lightColor = l.color();
	//PRINT_IN_MAYA(MString("Color is:") + colorToString(lightColor));

	LightDataT ld;
	ld.type = LightDataT::AMBIENT;
	ld.color = l.color();
	ld.intencity = l.intensity();
	lightingData.push_back(ld);
}

void RayTracer::storeDirectionalLight(MDagPath lightDagPath)
{
	MStatus status;
	MFnLight l(lightDagPath);

	MVector lightDir(0,0,-1); // origin
	lightDir *= lightDagPath.inclusiveMatrix();
	lightDir.normalize();

	LightDataT ld;
	ld.type = LightDataT::DIRECTIONAL;
	ld.color = l.color();
	ld.intencity = l.intensity();
	ld.direction = lightDir;
	lightingData.push_back(ld);
}

void RayTracer::storePointLight(MDagPath lightDagPath)
{
	MStatus status;
	MFnLight l(lightDagPath);

	LightDataT ld;
	ld.type = LightDataT::POINT;
	ld.color = l.color();
	ld.intencity = l.intensity();
	ld.position = MPoint(0,0,0,1) * lightDagPath.inclusiveMatrix();
	lightingData.push_back(ld);
}

#pragma endregion

void RayTracer::storeMeshMaterial(MeshDataT& m, const MDagPath& path)
{
	MFnMesh fn(path);
	MObjectArray shaders;
	MIntArray indices;
	fn.getConnectedShaders(0, shaders, indices);
	for(uint i = 0; i < shaders.length(); i++)
	{

		MPlugArray connections;
		MFnDependencyNode shaderGroup(shaders[i]);
		MPlug shaderPlug = shaderGroup.findPlug("surfaceShader");
		shaderPlug.connectedTo(connections, true, false);
		for(uint u = 0; u < connections.length(); u++)
		{
			if(connections[u].node().hasFn(MFn::kLambert))
			{
				MFnDependencyNode* dp = new MFnDependencyNode(connections[u].node());
				MStringArray sets;
				m.material.load(dp, sets);
				
			}
			else{
				//default material
				m.material.toDefault();
			}
		} 
	}
}

void RayTracer::computeAndStoreMeshData()
{
	MStatus status;
	MItDag dagIterator(MItDag::kDepthFirst, MFn::kMesh , &status);
	
	for(; !dagIterator.isDone(); dagIterator.next())
	{

		MDagPath dagPath;
		status = dagIterator.getPath(dagPath);


		triangulateMesh(MFnMesh(dagPath));

		pair<MPoint,MPoint> boundingBox = computeWfAxisAlignedBoundingBox(dagPath);
		
		MeshDataT aMesh;
		aMesh.min = boundingBox.first;
		aMesh.max = boundingBox.second;
		storeMeshMaterial(aMesh,dagPath);

		MFnMesh meshFn(dagPath);
		int faceCount = meshFn.numPolygons();

		totalPolyCount += faceCount; // Statistics

		vector<Face>& faces = aMesh.faces;
		faces.resize(faceCount);
		MItMeshPolygon faceIt(dagPath);
		bool isMeshTextured = aMesh.material.isTextured;
		for(int faceId = 0; !faceIt.isDone(); faceIt.next(), ++faceId)
		{
			Face& f = faces[faceId];
			faceIt.getPoints(f.vertices, MSpace::kWorld);
			faceIt.getNormals(f.normals, MSpace::kWorld);
			if(isMeshTextured)
			{
				faceIt.getUVs(f.us, f.vs);
			}
		}
		meshesData.push_back(aMesh); 

#ifdef PRINT_FOR_DEBUG
		PRINT_IN_MAYA(MString("Storing mesh, bb is:") + pointToString(aMesh.min) + "," + pointToString(aMesh.max));
#endif
	}

	MSelectionList selected;
	MGlobal::getActiveSelectionList(selected);
	selected.clear();
	MGlobal::setActiveSelectionList(selected);
}

void RayTracer::computeAndStoreSceneBoundingBox()
{
	minScene = MPoint( DBL_MAX ,DBL_MAX,DBL_MAX);
	maxScene = MPoint(-DBL_MAX, -DBL_MAX, -DBL_MAX);
	for (int i = 0; i < meshesData.size(); i++)
	{
		minimize(&(minScene.x), meshesData[i].min.x);
		minimize(&(minScene.y), meshesData[i].min.y);
		minimize(&(minScene.z), meshesData[i].min.z);

		maximize(&(maxScene.x), meshesData[i].max.x);
		maximize(&(maxScene.y), meshesData[i].max.y);
		maximize(&(maxScene.z), meshesData[i].max.z);
	}

	sceneBBPlanes.clear();
	sceneBBPlanes.resize(6);

	sceneBBPlanes[X_NEG] = Plane(minScene, MVector(-1,0,0));
	sceneBBPlanes[X_POS] = Plane(maxScene, MVector(1,0,0));
	sceneBBPlanes[Y_NEG] = Plane(minScene, MVector(0, -1 ,0));
	sceneBBPlanes[Y_POS] = Plane(maxScene, MVector(0,1,0));
	sceneBBPlanes[Z_NEG] = Plane(minScene, MVector(0,0, -1));
	sceneBBPlanes[Z_POS] = Plane(maxScene, MVector(0,0,1));

#ifdef PRINT_FOR_DEBUG
	PRINT_IN_MAYA(MString("Scene, bb is:") + pointToString(minScene) + "," + pointToString(maxScene));
#endif
}

void RayTracer::voxelizeScene()
{
	computeAndStoreVoxelParams();
	computeAndStoreRawVoxelsData();
	computeVoxelMeshIntersections();
}

void RayTracer::computeAndStoreVoxelParams()
{
	double sceneSpanX = maxScene.x - minScene.x;
	double sceneSpanY = maxScene.y - minScene.y;
	double sceneSpanZ = maxScene.z - minScene.z;

	sceneParams.dimensionDeltas[0] = sceneSpanX / sceneParams.voxelsPerDimension;
	sceneParams.dimensionDeltas[1] = sceneSpanY / sceneParams.voxelsPerDimension;
	sceneParams.dimensionDeltas[2] = sceneSpanZ / sceneParams.voxelsPerDimension;

	sceneParams.dimensionDeltaHalfs[0] = sceneParams.dimensionDeltas[0] / 2;
	sceneParams.dimensionDeltaHalfs[1] = sceneParams.dimensionDeltas[1] / 2;
	sceneParams.dimensionDeltaHalfs[2] = sceneParams.dimensionDeltas[2] / 2;
}

void RayTracer::computeAndStoreRawVoxelsData()
{
	for (int i = 0; i < voxelsData.size(); i++) {
		if (voxelsData[i].v != NULL) {
			delete voxelsData[i].v;
		}
	}
	voxelsData.clear();

	int sideCount = sceneParams.voxelsPerDimension;
	int maxIdx = sideCount * sideCount * sideCount - 1;
	voxelsData.resize(maxIdx + 1);


	double dx = sceneParams.dimensionDeltas[0];
	double dy = sceneParams.dimensionDeltas[1];
	double dz = sceneParams.dimensionDeltas[2];

	double x = minScene.x;
	double y = minScene.y;
	double z = minScene.z;

	cameraInSceneBB = isPointInVolume(activeCameraData.eye, minScene, maxScene);

	// TODO: Can reverse loops order to improve performance
	x = minScene.x;;
	for (int ix = 0; ix < sideCount; x += dx, ix++)
	{
		y = minScene.y;
		for (int iy = 0; iy < sideCount; y += dy, iy++)
		{
			z = minScene.z;
			for (int iz = 0; iz < sideCount; z += dz, iz++)
			{
				Voxel* v = new Voxel(MPoint(x,y,z), MPoint(x + dx, y + dy, z + dz) );
				VoxelDataT vd;
				vd.v = v;
				int index = sceneParams.flatten3dCubeIndex(ix, iy, iz);

				if(cameraInSceneBB && isPointInVolume(activeCameraData.eye, v->Min(), v->Max()))
				{
					initCameraVoxelX = ix;
					initCameraVoxelY = iy;
					initCameraVoxelZ = iz;
				}

				voxelsData[index] = vd;
			}
		}
	}
}

void RayTracer::computeVoxelMeshIntersections()
{
	double dimensionDeltaHalfs[3];
	for (int i = 0; i < 3; ++i)
	{
		dimensionDeltaHalfs[i] = sceneParams.dimensionDeltaHalfs[i] + DOUBLE_NUMERICAL_THRESHHOLD;
	}

	uint voxelNum = (uint)voxelsData.size();
	uint meshNum = (uint)meshesData.size();
	uint v, mid;
	for (v = 0; v < voxelNum; v++)
	{
		for (mid = 0; mid < meshNum; ++mid)
		{
			if ( voxelsData[v].v->intersectsWith(meshesData[mid].min, meshesData[mid].max) )
			{
				vector<int> faceIds;
				if( voxelsData[v].v->intersectsWith(meshesData[mid], dimensionDeltaHalfs, faceIds))
				{
					voxelsData[v].meshIdToFaceIds[mid] = faceIds;
				}
			}
		}
	}

}

void RayTracer::bresenhaim()
{
	int width = imagePlane.imgWidth;
	int height = imagePlane.imgHeight;
	int totalPixels = width*height;
	unsigned char* pixels = new unsigned char[totalPixels*4];
	double* pixelTimes = new double[totalPixels];
	int* pixelSamples =  new int[totalPixels];
	memset(pixels,0,totalPixels*4);
	memset(pixelTimes,0,totalPixels*sizeof(double));
	

#pragma region PARALLEL COMPUTATION
#pragma omp parallel for schedule(dynamic,100) num_threads(8)
	for(int it = 0; it < totalPixels; ++it)
	{
		MTimer timer;
		timer.beginTimer();
		int w = it % width;
		int h = it / width;
		MColor pixelColor;
		vector<MPoint> pointsOnPlane;
		
		imagePlane.getPointsOnIP(w, h, pointsOnPlane);

		switch (imagePlane.ssType) {
		case ImagePlaneDataT::UNIFORM:
		case ImagePlaneDataT::JITTERED:
		case ImagePlaneDataT::RANDOM:
			{
				int count = pointsOnPlane.size();
				MPoint raySource = activeCameraData.eye;
				MVector rayDirection = activeCameraData.viewDir;
				pixelSamples[it] = count;
				for(int ssit = 0; ssit < count; ++ssit )
				{
					if (activeCameraData.isPerspective) {
						rayDirection = (pointsOnPlane[ssit] - activeCameraData.eye).normal();
					}
					else {
						raySource = pointsOnPlane[ssit];
					}
					int depth = 0;
					pixelColor = sumColors(pixelColor, shootRay(raySource, rayDirection, sceneParams.rayDepth, &depth) / ((float)(count)));
#pragma omp atomic 
					totalDepths += depth;
				}
			}
			break;
		case ImagePlaneDataT::ADAPTIVE:
			bool needToStop = false;
			int count = 0;
			MPoint raySource = activeCameraData.eye;
			MVector rayDirection = activeCameraData.viewDir;

			MColor newColor;
			MColor colorExpectation;
			MColor colorPrevExpectation;
			MColor colorVariance;
			pixelSamples[it] = 0;
			while (!needToStop) {
				pixelSamples[it] ++;
				MPoint nextPoint = imagePlane.nextRandomPointOnIP(w, h);
				if (activeCameraData.isPerspective) {
					rayDirection = (nextPoint - activeCameraData.eye).normal();
				}
				else {
					raySource = nextPoint;
				}
				int depth = 0;
				newColor = shootRay(raySource, rayDirection, sceneParams.rayDepth, &depth); 
#pragma omp atomic 
				totalDepths += depth;
				count++; 
				
				if (1 == count) // first ray - here we initialize all the variance things
				{
					pixelColor = newColor;
					colorPrevExpectation = newColor;
					colorExpectation = newColor;
					colorVariance = MColor(0,0,0,0);
				}
				else 
				{
					pixelColor = nextColorAverage(pixelColor, count, newColor);
					colorPrevExpectation = colorExpectation;
					colorExpectation = nextColorExpectation(colorExpectation, count, newColor);
					colorVariance = nextColorVariance(colorVariance, colorPrevExpectation, colorExpectation, count, newColor);
				}

				if (count >= imagePlane.ssAdaptiveMaxSamples) 
				{
					needToStop = true;
				}
				else if (count >= imagePlane.ssAdaptiveMinSamples) {
					if (varianceIsSmallEnough(colorVariance, count, imagePlane.ssAdaptiveTolerance, imagePlane.ssAdaptiveErrorProbability)) {
						needToStop = true;
					}
				}
			}
			break;
		}

		pixels[h*width*4 + w*4] = (unsigned char) (pixelColor.r * 255.0);
		pixels[h*width*4 + w*4 + 1] = (unsigned char) (pixelColor.g * 255.0);
		pixels[h*width*4 + w*4 + 2] = (unsigned char) (pixelColor.b * 255.0);

		timer.endTimer();
		pixelTimes[it] = timer.elapsedTime();
	}
#pragma endregion

	computePixelStatistics(pixelTimes,pixelSamples, totalPixels);

	MImage img;
	img.setPixels(pixels,width,height);
	img.writeToFile(outputFilePath);
	img.release();
	delete [] pixels;
	delete [] pixelTimes;
}

void RayTracer::computePixelStatistics(double* pixelTimes,int* pixelSamples, int size)
{
	double sumTimePerPixel = 0;
	double averageTimePerPixel = 0;
	double varianceTimePerPixel = 0;

	double sumSamples = 0;
	double avgSamples = 0;
	double varSamples = 0;

	for (int i = 0; i < size; i++) 
	{
		sumTimePerPixel += pixelTimes[i];
		sumSamples += pixelSamples[i];
	}

	averageTimePerPixel = sumTimePerPixel / (double)size;
	avgSamples = sumSamples / (double) size;
	for (int i = 0; i < size; i++) 
	{
		varianceTimePerPixel += (pixelTimes[i] - averageTimePerPixel) * (pixelTimes[i] - averageTimePerPixel);
		varSamples += (pixelSamples[i] - avgSamples) * (pixelSamples[i] - avgSamples);
	}
	varianceTimePerPixel = varianceTimePerPixel/(double)size;
	varSamples = varSamples / (double)size;

	timePerPixel = averageTimePerPixel;
	timePerPixelStandardDeviation = sqrt(varianceTimePerPixel);

	samplesPerPixel = avgSamples;
	samplesPerPixelStdDeviation = sqrt(varSamples);
	totalSamples = sumSamples;
}

bool getOutRay(const MeshDataT& mesh, const MVector& view, const MPoint& inPoint, const MVector& inRay, MPoint& outPoint, MVector& outRay)
{
	MVector dir = inRay;
	MPoint src = inPoint + dir * DOUBLE_NUMERICAL_THRESHHOLD * 100;
	int size = mesh.faces.size();

	for( int count = 100; count > 0; --count) {
	
		bool intersected = false;
		MPoint tIntersection, mintIntersection;
		int outFaceId;
		double minTime = DBL_MAX, time = DBL_MAX;

		for(int fi = 0; fi < size; ++ fi) {
			if(rayIntersectsTriangle(src + dir * DOUBLE_NUMERICAL_THRESHHOLD * 100, dir, mesh.faces[fi].vertices, time, tIntersection) && time < minTime) {
				intersected = true;
				minTime = time;
				mintIntersection = tIntersection;
				outFaceId = fi;
			}
		}

		if( ! intersected) // out face is not found
			return false;

		const Face& outFace = mesh.faces[outFaceId];
		double bc2[3]; // baricentric coords
		calculateBaricentricCoordinates(outFace.vertices, mintIntersection, bc2 );

		MVector normal2 = (bc2[0] * outFace.normals[0] +  bc2[1] * outFace.normals[1] +  bc2[2] * outFace.normals[2]).normal();
		if(normal2* dir > 0 ) 
			normal2 = - normal2;
		MVector r;
		if( transmissionRay( dir, normal2, mesh.material.refractiveIndex, 1, r)) {
			outRay = r;
			outPoint = mintIntersection;
			return true;
		}	

		dir = reflectedRay(inRay, normal2);
		src = mintIntersection + dir * DOUBLE_NUMERICAL_THRESHHOLD * 100;

	}
	return false;
}


void RayTracer::calculateSpecularAndDiffuseCoeffs(const MPoint& intersection, const MVector& lightDir, const double distDepth, const MVector& normal, const MVector& view, int x, int y, int z, double& kd, double& ks) 
{ 
	kd = ks = 0.0;
	int meshId, faceId;
	MPoint secondIntersection;

	if(! closestIntersection(intersection, -lightDir , x, y, z, meshId , faceId, secondIntersection, distDepth )){
		kd = std::max(- (lightDir * normal), 0.0);
		ks = std::max( -(reflectedRay(lightDir, normal) * view) , 0.0);
	}
}

MColor RayTracer::shootRay(const MPoint& raySrc, const MVector& rayDir, int depth, int* depthReached)
{
	int x, y, z;
	int meshIdx, faceIdx;
	MPoint intersection;
	if (!findStartingVoxelIndeces(raySrc, rayDir, x, y, z) ||
		!closestIntersection(raySrc, rayDir, x, y, z, meshIdx, faceIdx, intersection )) {
			return BACKGROUND_COLOR;
	}

	MeshDataT& mesh = meshesData[meshIdx];
	Face& face = meshesData[meshIdx].faces[faceIdx];
	Material& mat = mesh.material;

#pragma region MeshPrecalculations
	double bc[3]; // baricentric coords
	calculateBaricentricCoordinates(face.vertices, intersection, bc );

	MVector normal = (bc[0] * face.normals[0] +  bc[1] * face.normals[1] +  bc[2] * face.normals[2]).normal();

	MColor diffuseMaterialColor;
	// TODO: diffuse coefficient
	if (!mat.isTextured) {
		diffuseMaterialColor = mat.diffuse;
	}
	else {
		// get texture color at point using u,v and bilinear filter
		double u = bc[0] * face.us[0] +  bc[1] * face.us[1] +  bc[2] * face.us[2];
		double v = bc[0] * face.vs[0] +  bc[1] * face.vs[1] +  bc[2] * face.vs[2];
		diffuseMaterialColor = getBilinearFilteredPixelColor(mat.texture, u, v);
	}


#pragma endregion

	MColor pixelColor = MColor(0,0,0,1);
	for (int li = (int) lightingData.size()- 1; li >= 0; --li)
	{
		LightDataT & currLight = lightingData[li];
		MColor lightColor = currLight.color * currLight.intencity;

		if( LightDataT::AMBIENT == currLight.type) {
			pixelColor = sumColors(mat.ambient * lightColor * (1 - mat.transparency), pixelColor);
		}
		else if(LightDataT::DIRECTIONAL == currLight.type || LightDataT::POINT == currLight.type)
		{
			double kd(0.0), ks(0.0);
			calculateSpecularAndDiffuseCoeffs(intersection, currLight.directionToPoint(intersection), currLight.distanceToPoint(intersection), normal, rayDir, 
				x, y, z, kd, ks);

			pixelColor = sumColors(pixelColor, diffuseMaterialColor * kd * lightColor * (1 - mat.transparency) * mat.diffuseCoeff);
			if( mat.cosPower > 1)
				pixelColor = sumColors(pixelColor, mat.specular * pow(ks, mat.cosPower) * lightColor);
		}
	}

	if (NULL != depthReached) {
		*depthReached = 1;
	}

	MVector inRay;
	if (depth < 1 || ! transmissionRay(rayDir , normal, 1, mesh.material.refractiveIndex, inRay))
		return pixelColor;

	double cos1 = - rayDir *  normal;
	double angle1 = acos(cos1);
	double cos2 = - inRay *  normal;
	double angle2 = acos(cos2);

	double sum = angle1 + angle2;
	double dif = angle1 - angle2;
	double kr = pow( sin(dif), 2)/ pow(sin(sum), 2) * ( 1 + pow(cos(sum), 2)/ pow( cos(dif), 2)) / 2;
	kr = std::min(1.0,kr);
	//double kr = mat.kr0 + (1 - mat.kr0) * pow( 1 - coss, 5);
	double kt = (1 - kr) * mat.transparency;
	
	double effectiveTransparency = ((1 - mat.reflectivity) * (1 - kr)) * mat.transparency;
	double effectiveReflectivity = mat.reflectivity + kr * (1 - mat.reflectivity);

	int transparentDepth = 0;
	int reflectedDepth = 0;

	if(mat.isTransparent){
		MVector outRay;
		MPoint outPoint;
		if(getOutRay(mesh, rayDir, intersection, inRay, outPoint, outRay)){
				MColor second = shootRay(outPoint, outRay, depth - 1, &transparentDepth);
				pixelColor = sumColors( pixelColor , second * effectiveTransparency);
		}
	}
	if(mat.isReflective) {
		MVector reflected = reflectedRay(rayDir, normal);
		MColor reflColor = shootRay(intersection + reflected * 0.001, reflected, depth - 1, &reflectedDepth);
		pixelColor = sumColors(pixelColor, reflColor * effectiveReflectivity ); 
	}
	if (NULL != depthReached) {
		*depthReached = 1 + std::max(transparentDepth, reflectedDepth);
		if (*depthReached > 2) {
			bool ok = true;
		}
	}

	return pixelColor;
}

bool RayTracer::findStartingVoxelIndeces(const MPoint& raySrc, const MVector& rayDirection, int& x, int& y, int& z)
{

	//if(activeCameraData.isPerspective && cameraInSceneBB)
	//{
	//	x = initCameraVoxelX;
	//	y = initCameraVoxelY;
	//	z = initCameraVoxelZ;
	//	return true;
	//}
	x = y = z = 0;
	if(isPointInVolume(raySrc, minScene, maxScene)  &&
		findIndecesByDimension( raySrc, X_POS, x, y, z) &&
		findIndecesByDimension( raySrc, Y_POS, x, y, z) && 
		findIndecesByDimension( raySrc, Z_POS, x, y, z))
	{
		return true;
	}

	MPoint closestIntersection;
	double time = DBL_MAX;
	AxisDirection direction = UNKNOWN_DIR;

	for (int i = 0; i < UNKNOWN_DIR; ++i)
	{
		MPoint curIntersection;
		double curTime;
		AxisDirection curDirection = (AxisDirection)i;
		if( sceneBBPlanes[i].rayIntersection(raySrc, rayDirection, curTime, curIntersection)
			&& pointInRectangle(curDirection, curIntersection, minScene, maxScene))
		{
			if(curTime < time)
			{
				closestIntersection = curIntersection;
				time = curTime;
				direction = curDirection;
			}
		}
	}

	if(direction == UNKNOWN_DIR)
		return false;

	initIndeces(direction, x,y,z);
	AxisDirection uDirection, vDirection;
	orthonormalDirections(direction, uDirection, vDirection);

	if (!findIndecesByDimension( closestIntersection, uDirection, x, y, z) ||
		!findIndecesByDimension(closestIntersection, vDirection, x, y, z) )
	{
		return false;
	}
	return true;
}

bool RayTracer::findIndecesByDimension( const MPoint& point, AxisDirection direction, int& x, int& y, int& z )
{
	int dimension = sceneParams.voxelsPerDimension;
	int cur3DIndex = sceneParams.flatten3dCubeIndex(x,y,z);
	for(;x >= 0 && x < dimension && y >= 0 && y < dimension && z >= 0 && z < dimension; sceneParams.incrementIndeces(direction, x, y, z, cur3DIndex))
	{
		if(pointInVoxelByDirection(point, voxelsData[cur3DIndex], direction))
		{
			return true;
		}
	}
	return false;
}

void RayTracer::initIndeces( AxisDirection direction, int& x, int& y, int& z )
{
	x = 0;
	y = 0;
	z = 0;
	if (direction == X_POS)
	{
		x = sceneParams.voxelsPerDimension - 1;
	}
	else if (direction == Y_POS)
	{
		y = sceneParams.voxelsPerDimension - 1;

	}
	else if (direction == Z_POS)
	{
		z = sceneParams.voxelsPerDimension - 1;
	}
}

void RayTracer::orthonormalDirections( AxisDirection direction, AxisDirection& uDirection, AxisDirection& vDirection )
{
	switch (direction)
	{
	case X_NEG:
	case X_POS:
		uDirection = Y_POS;
		vDirection = Z_POS;
		break;
	case Y_NEG:
	case Y_POS:
		uDirection = X_POS;
		vDirection = Z_POS;
		break;
	case Z_NEG:
	case Z_POS:
		uDirection = X_POS;
		vDirection = Y_POS;
		break;
	}
}

inline bool RayTracer::pointInVoxelByDirection( const MPoint& closestIntersection, VoxelDataT& voxel, AxisDirection direction )
{
	switch (direction)
	{
	case X_POS:
		return valueInInterval(closestIntersection.x, voxel.v->Min().x, voxel.v->Max().x);
	case Y_POS:
		return valueInInterval(closestIntersection.y, voxel.v->Min().y, voxel.v->Max().y);
	case Z_POS:
		return valueInInterval(closestIntersection.z, voxel.v->Min().z, voxel.v->Max().z);
	default:
		return false;
	}

}

// The function finds the mesh with it intersects given ray, the inner id of the face in the mesh and the intersection point.
// Also it changes the x,y,z indices to match the voxel where the closest intersection happens.
// Returns true if finds
// Return false if it arrives to the scene bounds and doesn't meet any mesh an some point.
bool RayTracer::closestIntersection(const MPoint& raySource,const MVector& rayDirection, int& x, int& y, int& z , int& meshIndex, int& innerFaceId, MPoint& intersection , double depth)
{
#pragma omp atomic
	totalRayCount++;

	AxisDirection  farAxisDir;
	int currMeshIndex, currInnerFaceId;
	MPoint currIntersection;
	int cur3Dindex = sceneParams.flatten3dCubeIndex( x, y, z);

	for(	;
			x >= 0 && x < sceneParams.voxelsPerDimension && y >= 0 && y < sceneParams.voxelsPerDimension && z >= 0 && z < sceneParams.voxelsPerDimension; 
			sceneParams.incrementIndeces(farAxisDir, x, y, z, cur3Dindex)) 
	{
		VoxelDataT& voxelData = voxelsData[cur3Dindex];
		if(!voxelData.v->findExitDirection(raySource, rayDirection, farAxisDir)) {
			break;
		}
		if(voxelData.meshIdToFaceIds.size() == 0) {
			continue;
		}
		if(!closestIntersectionInVoxel(raySource, rayDirection, voxelData, currMeshIndex, currInnerFaceId, currIntersection)) {
			continue;
		}
		double distToIntersection = (currIntersection - raySource).length();
		if( distToIntersection > depth ) {
			return false;
		}
		if (distToIntersection < DOUBLE_NUMERICAL_THRESHHOLD) {
			continue;
		}
		meshIndex = currMeshIndex;
		innerFaceId = currInnerFaceId;
		intersection = currIntersection;
		return true;
	}

	return false;
}

bool RayTracer::closestIntersectionInVoxel(const MPoint& raySource, const MVector& rayDirection, VoxelDataT &voxelData, int &meshIndex, int &innerFaceId, MPoint &intersection )
{

	bool res = false;
	double minTime = DBL_MAX;
	double time;
	int currentFaceIndex;
	MPoint curIntersection;

	for(map<int, vector<int>>::iterator it = voxelData.meshIdToFaceIds.begin(); it != voxelData.meshIdToFaceIds.end(); ++it )
	{
		MeshDataT& mesh = meshesData[it->first];
		vector<int>& faceIds = it->second;
		for(currentFaceIndex = (int) faceIds.size() - 1; currentFaceIndex >= 0; --currentFaceIndex)
		{
#pragma omp atomic
			intersectionTestCount++;
			
			Face& face = mesh.faces[faceIds[currentFaceIndex]];
			/*
			mesh.getPolygonVertices(faceIds[currentFaceIndex], vertexIds);
			if(vertexIds.length() != 3) {
			continue;
			}
			for (int vi = 0; vi < 3; ++vi) {
			mesh.getPoint(vertexIds[vi], triangleVertices[vi], MSpace::kWorld);
			}*/

			if(!rayIntersectsTriangle(raySource, rayDirection, face.vertices, time, curIntersection)
				|| !isPointInVolume(curIntersection, voxelData.v->Min(), voxelData.v->Max())
				|| ((curIntersection - raySource)*rayDirection) < 0)
			{
				continue;
			}
			if(time < minTime) {
				meshIndex = it->first;
				innerFaceId = faceIds[currentFaceIndex];
				intersection = curIntersection;
				minTime = time;
				res = true;
#pragma omp atomic
			intersectionFoundCount++;
				
			}
		}
	}

	return res;
}



#include "RayTracer.h"

#include "Profiler.h"


#ifdef _DEBUG
#define DEBUG_REPORT 0
#endif

#ifdef PRINT_FOR_DEBUG
#define PRINT_IN_MAYA(arg) MGlobal::displayInfo((arg));
#endif

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


char*	RayTracer::outputFilePath = "C://temp//scene.iff";
char*	RayTracer::statisticsFilePath = "C://temp//stat.txt";


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

	syntax.addFlag(widthFlag, "-width", MSyntax::kLong);
	syntax.addFlag(heightFlag, "-height", MSyntax::kLong);
	syntax.addFlag(voxelsFlag, "-voxels", MSyntax::kLong);
	syntax.addFlag(supersamplingFlag, "-supersampling", MSyntax::kLong);

	return syntax;
}

bool RayTracer::parseArgs( const MArgList& args)
{
	MArgParser    argData(syntax(), args);
	MStatus s;
	MString         arg;

	if ( argData.isFlagSet(widthFlag) ) {
		uint arg;
		s = argData.getFlagArgument(widthFlag, 0, arg);	
		if (s == MStatus::kSuccess) {
			imgWidth = (arg < 1) ? 1 : arg;
		}
	}

	if ( argData.isFlagSet(heightFlag) ) {
		uint arg;
		s = argData.getFlagArgument(heightFlag, 0, arg);	
		if (s == MStatus::kSuccess) {
			imgHeight = (arg < 1) ? 1 : arg;
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
			// do stuff
			supersamplingCoeff = (arg < 1) ? 1 : arg;
		}
	}

	return true;
}

RayTracer::RayTracer()
{
	imgWidth = 1920;
	imgHeight = 1080;
	sceneParams.voxelsPerDimension = 1;
	supersamplingCoeff = 1;
	sceneParams.voxelsPerDimensionSqr = 1;
	
	
	prepTime = 0;
	totalTime = 0;
	timePerPixel = 0;
	timePerPixelStandardDeviation = 0;
	intersectionTestCount = 0;
	intersectionFoundCount = 0;
	voxelsTraversed = 0;
	totalRayCount = 0;
	totalPolyCount = 0;

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
	activeCameraData.focalLengthCm	= (camera.focalLength() / 10);

	if (camera.isOrtho() ) {
		activeCameraData.isPerspective	= false;
		activeCameraData.filmWidthCm	= camera.orthoWidth();
	}
	else {
		activeCameraData.isPerspective	= true;
		activeCameraData.filmWidthCm	= 2.54 * camera.horizontalFilmAperture();
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
	MVector xAxis = activeCameraData.viewDir ^ activeCameraData.upDir;
	xAxis.normalize();
	MVector yAxis = xAxis ^ activeCameraData.viewDir;
	yAxis.normalize();

	imagePlane.x = xAxis;
	imagePlane.y = yAxis;

	MPoint centerPoint = activeCameraData.eye + ( activeCameraData.viewDir * activeCameraData.focalLengthCm );


	double imgAspect = (double)imgWidth / (double)imgHeight;
	//double pixelWidth = (activeCameraData.filmWidthCm)/(double)imgWidth;
	//double pixelHeight = pixelWidth / imgAspect;
	imagePlane.dp = ( activeCameraData.filmWidthCm)/(double)imgWidth;

	int halfWidth = imgWidth / 2;
	int halfHeight = imgHeight / 2;

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

void RayTracer::storeMeshTexturingData(MeshDataT& m)
{
	MFnMesh fn(m.dagPath);
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
				MFnLambertShader lambertShader(connections[u].node());
				MImage* textureImage = new MImage();
				if (getLambertShaderTexture(lambertShader, *textureImage)) 
				{
					// store the texture image
					m.texture = textureImage;
					m.hasTexture = true;
				}
				else 
				{
					// assign diffuse and ambient color
					m.hasTexture = false;
					m.diffuse = lambertShader.color();
				}
				m.useHalfVector = false;
				m.ambient = lambertShader.ambientColor();
				m.specularPower = -1.0;
				if(connections[u].node().hasFn(MFn::kReflect))
				{
					MFnReflectShader shader(connections[u].node());
					m.specular = shader.specularColor();
				}


				if(connections[u].node().hasFn(MFn::kPhong)) {
					MFnPhongShader shader(connections[u].node());
					m.specularPower = shader.cosPower();
				}
				else if(connections[u].node().hasFn(MFn::kBlinn)) {
					MFnBlinnShader blinn(connections[u].node());
					m.specularPower = 20;
					m.specular = sumColors(MColor(0,0,0), m.specular * blinn.specularRollOff() );
					m.useHalfVector = true;
					m.eccentricity = blinn.eccentricity();
				}
			}
			else{
				//default material
				m.ambient = MColor(0.1f, 0.1f, 0.1f);
				m.diffuse = MColor(0.f, 0.f, 0.8f);
				m.specular = MColor(0.9f,0.9f,0.9f);
				m.specularPower = 10;
				m.useHalfVector = false;
				m.hasTexture = false;
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
		CHECK_MSTATUS(status);
		
		triangulateMesh(MFnMesh(dagPath));
		pair<MPoint,MPoint> boundingBox = computeWfAxisAlignedBoundingBox(dagPath);
		MeshDataT aMesh;

		aMesh.dagPath = dagPath;
		aMesh.min = boundingBox.first;
		aMesh.max = boundingBox.second;
		storeMeshTexturingData(aMesh);

		MItMeshPolygon faceIt(dagPath);
		for(; !faceIt.isDone(); faceIt.next())
		{
			Face f;
			faceIt.getPoints(f.vertices, MSpace::kWorld);
			faceIt.getNormals(f.normals, MSpace::kWorld);
			if(aMesh.hasTexture)
			{
				faceIt.getUVs(f.us, f.vs);
			}
			aMesh.faces.push_back(f);
		}

		meshesData.push_back(aMesh);

		MFnMesh meshFn(dagPath);

		totalPolyCount += meshFn.numPolygons();

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

	int totalIntersections = 0;
	for (int v = 0; v < voxelsData.size(); v++)
	{
		int count = 0;
		voxelsData[v].containedMeshIndexes.clear();
		for (int m = 0; m < meshesData.size(); m++)
		{
			if ( voxelsData[v].v->intersectsWith(meshesData[m].min, meshesData[m].max) )
			{
				vector<int> faceIds;
				if( voxelsData[v].v->intersectsWith(meshesData[m].dagPath, dimensionDeltaHalfs, faceIds))
				{
					voxelsData[v].meshIdToFaceIds[m] = faceIds;
				}
				voxelsData[v].containedMeshIndexes.push_back(m);
				count++;
				totalIntersections++;
			}
		}
		/*MString out = "Intersections per voxel ";
		out += pointToString(voxelsData[v].v->min);
		out += "x";
		out += pointToString(voxelsData[v].v->max);
		out += " is ";
		out += count;
		PRINT_IN_MAYA( out );*/
	}
#ifdef PRINT_FOR_DEBUG
	PRINT_IN_MAYA((MString("TOTAL: ") + totalIntersections));
#endif
}

void RayTracer::bresenhaim()
{
	unsigned char* pixels = new unsigned char[imgWidth*imgHeight*4];
	memset(pixels,0,imgWidth*imgHeight*4);


	double ssDp = imagePlane.dp / (double) (supersamplingCoeff + 1);
	MVector ssdx = imagePlane.x * ssDp;
	MVector ssdy = imagePlane.y * ssDp;

	MVector dx = imagePlane.x * imagePlane.dp;
	MVector dy = imagePlane.y * imagePlane.dp;
	

	MPoint leftBottomOfPixel,bottomInPixel, leftBottomInPixel;
	MPoint bottomOfPixel = imagePlane.lb;
	int x, y, z;
	
	int meshIntersectionIndex, innerFaceIntersectionIndex;
	MPoint intersectionPoint;

	MPoint raySource = activeCameraData.eye;
	MVector rayDirection = activeCameraData.viewDir;
	MVector orthoFilOffset = (activeCameraData.focalLengthCm * activeCameraData.viewDir);

	int pixelCountSoFar = 0;
	double timePerPixelMean = 0;
	double M2 = 0;
	int totalPixels = imgHeight * imgWidth;
	for( int h = 0; h < imgHeight; ++h,bottomOfPixel += dy )
	{
		leftBottomOfPixel = bottomOfPixel;
		for( int w = 0; w < imgWidth; ++w , leftBottomOfPixel += dx)
		{
			Profiler::startTimer("bresenhaim::timePerPixel");
			
			MColor pixelColor;
			
			bottomInPixel = leftBottomOfPixel + ssdy + ssdx;
			for (int sh = 1; sh <= supersamplingCoeff; sh++, bottomInPixel += ssdy )
			{
				leftBottomInPixel = bottomInPixel;
				for (int sw = 1; sw <= supersamplingCoeff; sw++, leftBottomInPixel += ssdx)
				{
					if (activeCameraData.isPerspective) {
						rayDirection = (leftBottomInPixel - activeCameraData.eye).normal();
					}
					else {
						raySource = leftBottomInPixel - orthoFilOffset;
					}

					bool foundStartingVoxel = findStartingVoxelIndeces(raySource, rayDirection, x, y, z);

					if (!foundStartingVoxel) {
						continue;
					}

					bool foundIntersection = closestIntersection(raySource, rayDirection, x, y, z, meshIntersectionIndex, innerFaceIntersectionIndex, intersectionPoint );

					if(!foundIntersection) {
						// Put the background
						// currently black
						continue;
					}

					pixelColor = sumColors(pixelColor, (calculatePixelColor(x, y, z, rayDirection, meshIntersectionIndex, innerFaceIntersectionIndex, intersectionPoint) / (float)(supersamplingCoeff*supersamplingCoeff)));

				}
			}
			pixels[h*imgWidth*4 + w*4] = (unsigned char) (pixelColor.r * 255.0);
			pixels[h*imgWidth*4 + w*4 + 1] = (unsigned char) (pixelColor.g * 255.0);
			pixels[h*imgWidth*4 + w*4 + 2] = (unsigned char) (pixelColor.b * 255.0);

			pixelCountSoFar++;
			double elapsed = Profiler::finishTimer("bresenhaim::timePerPixel");
			double delta = elapsed - timePerPixelMean;
			timePerPixel += elapsed;
			timePerPixelMean = timePerPixelMean + delta/(double)pixelCountSoFar;
			M2 += delta * (elapsed - timePerPixelMean);

			for (int cni = 1; cni < 10; cni++ ) {
				if (cni*totalPixels/10 == pixelCountSoFar) {
					cout << cni*10 << "% done" << endl;
				}
			}
			
		}
	}
	timePerPixel = timePerPixel / ( (double)pixelCountSoFar );
	timePerPixelStandardDeviation = sqrt( M2 / ((double)pixelCountSoFar - 1) );
	MImage img;
	img.setPixels(pixels,imgWidth,imgHeight);
	img.writeToFile(outputFilePath);
	img.release();
	delete [] pixels;
}

bool RayTracer::findStartingVoxelIndeces(const MPoint& raySrc, const MVector& rayDirection, int& x, int& y, int& z)
{
	
	if(activeCameraData.isPerspective && cameraInSceneBB)
	{
		x = initCameraVoxelX;
		y = initCameraVoxelY;
		z = initCameraVoxelZ;
		return true;
	}
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
bool RayTracer::closestIntersection(const MPoint& raySource,const MVector& rayDirection, int& x, int& y, int& z , int& meshIndex, int& innerFaceId, MPoint& intersection )
{
	totalRayCount++;

	AxisDirection  farAxisDir;
	int currMeshIndex, currInnerFaceId;
	MPoint currIntersection;
	int cur3Dindex = sceneParams.flatten3dCubeIndex( x, y, z);
	for(;x >= 0 && x < sceneParams.voxelsPerDimension &&
		 y >= 0 && y < sceneParams.voxelsPerDimension && 
		 z >= 0 && z < sceneParams.voxelsPerDimension; 
				sceneParams.incrementIndeces(farAxisDir, x, y, z, cur3Dindex)) {
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
				|| !isPointInVolume(curIntersection, voxelData.v->Min(), voxelData.v->Max()))
			{
				continue;
			}
			if(time < minTime) {
				meshIndex = it->first;
				innerFaceId = faceIds[currentFaceIndex];
				intersection = curIntersection;
				minTime = time;
				res = true;
				intersectionFoundCount++;
			}
		}
	}

	return res;
}

MColor RayTracer::calculatePixelColor(const int x, const int y, const int z,const MVector& rayDirection, const int meshIndex,const int innerFaceId,const MPoint& intersection)
{
	MeshDataT& mesh = meshesData[meshIndex];
	Face& face = meshesData[meshIndex].faces[innerFaceId];

	/*MIntArray vertexIds;
	MPoint triangleVertices[3];
	MVector triangleNormals[3];
	float us[3];
	float vs[3];
	mesh.getPolygonVertices(innerFaceId, vertexIds);
	if(vertexIds.length() != 3) {
		MColor material( 0.25,0.25,0.25 );
		return material;
	}
	for (int vi = 0; vi < 3; ++vi) {
		mesh.getPoint(vertexIds[vi], triangleVertices[vi], MSpace::kWorld);
		mesh.getFaceVertexNormal(innerFaceId, vertexIds[vi], triangleNormals[vi], MSpace::kWorld );
		triangleNormals[vi].normalize();
		mesh.getPolygonUV(innerFaceId, vi, us[vi], vs[vi]);
	}*/


	double bc[3]; // baricentric coords

	caclulateBaricentricCoordinates(face.vertices, intersection, bc );


	MColor	specularMaterialColor	=	mesh.specular;
	MColor	ambientMaterialColor	=	mesh.ambient;
	float	specularPower			=	mesh.specularPower;
	MColor	diffuseMaterialColor;
	bool useHalfVector = mesh.useHalfVector;
	float eccentricity = mesh.eccentricity;
	if (!mesh.hasTexture) {
		diffuseMaterialColor = mesh.diffuse;
	}
	else {
		// get texture color at point using u,v and bilinear filter
		double u = bc[0] * face.us[0] +  bc[1] * face.us[1] +  bc[2] * face.us[2];
		double v = bc[0] * face.vs[0] +  bc[1] * face.vs[1] +  bc[2] * face.vs[2];
		diffuseMaterialColor = getBilinearFilteredPixelColor(mesh.texture, u, v); 
	}
	MVector normalAtPoint = (bc[0] * face.normals[0] +  bc[1] * face.normals[1] +  bc[2] * face.normals[2]).normal();

	MColor pixelColor = MColor(0,0,0,1);

	int currX, currY, currZ;
	MPoint secondIntersection;
	int secondIntersectionMeshIndex, secondIntersectionFaceId;
	for (int li = (int) lightingData.size()- 1; li >= 0; --li)
	{
		LightDataT & currLight = lightingData[li];
		MColor mixedDiffuse = diffuseMaterialColor * currLight.color * currLight.intencity;
		MColor mixedSpecular = specularMaterialColor * currLight.color * currLight.intencity;	
		currX = x;
		currY = y;
		currZ = z;

		switch (currLight.type)
		{
		case LightDataT::AMBIENT:
			pixelColor = sumColors(ambientMaterialColor * currLight.color * currLight.intencity, pixelColor);
			break;
		case LightDataT::DIRECTIONAL:
			if(! closestIntersection(intersection, -currLight.direction, currX, currY, currZ, secondIntersectionMeshIndex, secondIntersectionFaceId, secondIntersection )){
				pixelColor = sumColors(calculateSpecularAndDiffuse(rayDirection, currLight.direction, normalAtPoint, mixedDiffuse, mixedSpecular, specularPower, useHalfVector, eccentricity), pixelColor);
			}
			break;
		case LightDataT::POINT:
			{
				MVector lightDirection = intersection - currLight.position;
				MVector lightDirectionNormalized = lightDirection.normal();

				if(closestIntersection(intersection, -lightDirectionNormalized, currX, currY, currZ, secondIntersectionMeshIndex, secondIntersectionFaceId, secondIntersection )){
					if((secondIntersection - intersection).length() >= lightDirection.length()) // intersection after the light
					{
						pixelColor = sumColors(calculateSpecularAndDiffuse(rayDirection, lightDirectionNormalized, normalAtPoint, mixedDiffuse, mixedSpecular, specularPower, useHalfVector, eccentricity), pixelColor) ;
					}
				}
				else 
				{
					pixelColor = sumColors(calculateSpecularAndDiffuse(rayDirection, lightDirectionNormalized, normalAtPoint, mixedDiffuse, mixedSpecular, specularPower, useHalfVector, eccentricity), pixelColor);
				}
			}
			break;
		default:
			break;
		}
	}

	return pixelColor;
}

MColor RayTracer::calculateSpecularAndDiffuse(const MVector& viewDirection, MVector& lightDirection,  MVector& normalAtPoint, MColor& mixedDiffuse, MColor& mixedSpecular, float specularPower, bool useHalfVector, float eccentricity)
{
	// diffuse
	MColor currColorComponent;
	float k = (float) (lightDirection * normalAtPoint);
	if (k < 0) {
		currColorComponent = sumColors(-k * mixedDiffuse, currColorComponent);
	}

	// specular
	if (specularPower > 0)
	{
		if (useHalfVector) {
			k = (float) ( (halfVector( -lightDirection, -viewDirection ) * normalAtPoint));
		}
		else {
			k = -(float)(reflectedRay(lightDirection, normalAtPoint) * viewDirection);
		}
		if(k > 0) {
			currColorComponent = sumColors(pow(k, specularPower) * mixedSpecular, currColorComponent);
		}	
	}
	
	return currColorComponent;
}

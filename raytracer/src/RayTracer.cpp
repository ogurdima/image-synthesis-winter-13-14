#include "RayTracer.h"

#include "Profiler.h"

#define PRINT_IN_MAYA(arg) MGlobal::displayInfo((arg));


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
	Profiler::finishTimer("doIt::prepTime");

	bresenhaim();

	Profiler::finishTimer("doIt::totalTime");
	Profiler::printReport();


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
			voxelParams.voxelsPerDimension = (arg < 1) ? 1 : arg;
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
	imgWidth = 640;
	imgHeight = 480;
	voxelParams.voxelsPerDimension = 10;
	supersamplingCoeff = 1;
	outputFilePath = "C://temp//scene.iff";

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

bool badMStatus(const MStatus& status, const MString& error)
{
	if(status != MStatus::kSuccess)
	{
		status.perror(error);
		return true;
	}
	return false;
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

void RayTracer::storeActiveCameraData()
{
	MDagPath cameraPath;
	MStatus status;
	M3dView::active3dView().getCamera( cameraPath );
	MFnCamera camera(cameraPath, &status);
	CHECK_MSTATUS(status);
	MVector up = camera.upDirection(MSpace::kWorld);
	up.normalize();
	MVector view = camera.viewDirection(MSpace::kWorld);
	view.normalize();
	MVector eye = camera.eyePoint(MSpace::kWorld, &status);
	CHECK_MSTATUS(status);
	double focalMm = camera.focalLength();
	double horizontalAperture = camera.horizontalFilmAperture();

	activeCameraData.eye = eye;
	activeCameraData.filmWidthCm = horizontalAperture;
	activeCameraData.focalLengthCm = (focalMm / 10);
	activeCameraData.upDir = up;
	activeCameraData.viewDir = view;

	PRINT_IN_MAYA(activeCameraData.toString());
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
	imagePlane.dp = (2.54 * activeCameraData.filmWidthCm)/(double)imgWidth;

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
			PRINT_IN_MAYA("Unsupported light");
		}
	}

	for (int i = 0; i < lightingData.size(); i++)
	{
		PRINT_IN_MAYA(lightingData[i].toString());
	}
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
	MFloatVector dirFloat = l.lightDirection(lightDagPath.instanceNumber(), MSpace::kWorld, &status);
	CHECK_MSTATUS(status);

	MVector vecc(0,0,1);

	MVector lightDir(0,0,-1); // origin
	lightDir *= lightDagPath.inclusiveMatrix();
	lightDir.normalize();
	MVector dir;
	dir.x = (double) dirFloat.x;
	dir.y = (double) dirFloat.y;
	dir.z = (double) dirFloat.z;
	dir.normalize();
	//PRINT_IN_MAYA(MString("Directional Light in WF:") + vectorToString(dir));

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
	MMatrix transform = lightDagPath.inclusiveMatrix(); /*getDagPathTransformationMatrix(lightDagPath, &status);*/
	CHECK_MSTATUS(status);
	MPoint pointLightOrigin = MPoint(0,0,0,1);
	MPoint positionWorldFrame = pointLightOrigin * transform;


	//PRINT_IN_MAYA(MString("Point Light in WF:") + pointToString(positionWorldFrame));

	LightDataT ld;
	ld.type = LightDataT::POINT;
	ld.color = l.color();
	ld.intencity = l.intensity();
	ld.position = positionWorldFrame;
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
				m.ambient = lambertShader.ambientColor();

				if(connections[u].node().hasFn(MFn::kReflect))
				{
					MFnReflectShader shader(connections[u].node());
					m.specular = shader.specularColor();
				}

				if(connections[u].node().hasFn(MFn::kPhong))
				{
					MFnPhongShader shader(connections[u].node());
					m.specularPower = shader.cosPower();
				}
				else {
					m.specularPower = -1.0;
				}

			}
		}
	}
}

void RayTracer::computeAndStoreMeshData()
{
	MStatus status;
	MItDag dagIterator(MItDag::kDepthFirst, MFn::kMesh , &status);
	CHECK_MSTATUS(status);
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
		meshesData.push_back(aMesh);

		MFnMesh meshFn(dagPath);
		Profiler::increaseCounter("totalPolygons", meshFn.numPolygons());
		PRINT_IN_MAYA(MString("Storing mesh, bb is:") + pointToString(aMesh.min) + "," + pointToString(aMesh.max));
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

	PRINT_IN_MAYA(MString("Scene, bb is:") + pointToString(minScene) + "," + pointToString(maxScene));
}

void RayTracer::voxelizeScene()
{
	computeAndStoreVoxelParams();
	computeAndStoreRawVoxelsData();
	computeVoxelMeshBboxIntersections();
}

void RayTracer::computeAndStoreVoxelParams()
{
	double sceneSpanX = maxScene.x - minScene.x;
	double sceneSpanY = maxScene.y - minScene.y;
	double sceneSpanZ = maxScene.z - minScene.z;

	//maximize(&sceneSpanX, 1); // to avoid zero voxel size in case of a flat scene
	//maximize(&sceneSpanY, 1); // to avoid zero voxel size in case of a flat scene
	//maximize(&sceneSpanZ, 1); // to avoid zero voxel size in case of a flat scene

	voxelParams.dimensionDeltas[0] = sceneSpanX / voxelParams.voxelsPerDimension;
	voxelParams.dimensionDeltas[1] = sceneSpanY / voxelParams.voxelsPerDimension;
	voxelParams.dimensionDeltas[2] = sceneSpanZ / voxelParams.voxelsPerDimension;

	voxelParams.dimensionDeltaHalfs[0] = voxelParams.dimensionDeltas[0] / 2;
	voxelParams.dimensionDeltaHalfs[1] = voxelParams.dimensionDeltas[1] / 2;
	voxelParams.dimensionDeltaHalfs[2] = voxelParams.dimensionDeltas[2] / 2;

	//PRINT_IN_MAYA((MString("Voxel size: ") + pointToString(MPoint(voxelParams.dx, voxelParams.dy, voxelParams.dz))));

}

void RayTracer::computeAndStoreRawVoxelsData()
{
	for (int i = 0; i < voxelsData.size(); i++) {
		if (voxelsData[i].v != NULL) {
			delete voxelsData[i].v;
		}
	}
	voxelsData.clear();

	int sideCount = voxelParams.voxelsPerDimension;
	int maxIdx = sideCount * sideCount * sideCount - 1;
	voxelsData.resize(maxIdx + 1);


	double dx = voxelParams.dimensionDeltas[0];
	double dy = voxelParams.dimensionDeltas[1];
	double dz = voxelParams.dimensionDeltas[2];

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
				int index = flatten3dCubeIndex(sideCount, ix, iy, iz);

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

void RayTracer::computeVoxelMeshBboxIntersections()
{
	double dimensionDeltaHalfs[3];
	for (int i = 0; i < 3; ++i)
	{
		dimensionDeltaHalfs[i] = voxelParams.dimensionDeltaHalfs[i] + DOUBLE_NUMERICAL_THRESHHOLD;
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
	PRINT_IN_MAYA((MString("TOTAL: ") + totalIntersections));
}

void RayTracer::bresenhaim()
{
	unsigned char* pixels = new unsigned char[imgWidth*imgHeight*4];
	memset(pixels,0,imgWidth*imgHeight*4);
	MPoint lbPixelCenter = imagePlane.lb + (imagePlane.x + imagePlane.y) * (imagePlane.dp / 2);

	double ssDp = imagePlane.dp / (double) (supersamplingCoeff + 1);
	MVector ssdx = imagePlane.x * ssDp;
	MVector ssdy = imagePlane.y * ssDp;

	MVector dx = imagePlane.x * imagePlane.dp;
	MVector dy = imagePlane.y * imagePlane.dp;

	int dimension = voxelParams.voxelsPerDimension;
	MPoint raySource = activeCameraData.eye;

	
	for( int h = 0; h < imgHeight; ++h )
	{
		for( int w = 0; w < imgWidth; ++w )
		{
			Profiler::startTimer("bresenhaim::timePerPixel");
			int meshIntersectionIndex, innerFaceIntersectionIndex;
			MPoint intersectionPoint;
			int x, y, z;
			MColor pixelColor;
			MPoint rayOnImagePlane;
			MVector rayDirection;

			for (int sh = 1; sh <= supersamplingCoeff; sh++ )
			{
				for (int sw = 1; sw <= supersamplingCoeff; sw++)
				{
					rayOnImagePlane = imagePlane.lb + ( h * dy ) + (w * dx) + (sh * ssdy) + (sw * ssdx);
					rayDirection = (rayOnImagePlane - activeCameraData.eye).normal();


					Profiler::startTimer("bresenhaim::findStartingVoxel");
					bool foundStartingVoxel = findStartingVoxelIndeces(rayDirection, x, y, z);
					Profiler::finishTimer("bresenhaim::findStartingVoxel");
					if (!foundStartingVoxel) {
						continue;
					}

					Profiler::startTimer("bresenhaim::closestIntersection");
					bool foundIntersection = closestIntersection(dimension, raySource, rayDirection, x, y, z, meshIntersectionIndex, innerFaceIntersectionIndex, intersectionPoint );
					Profiler::finishTimer("bresenhaim::closestIntersection");
					if(!foundIntersection) {
						// Put the background
						// currently black
						continue;
					}
					Profiler::startTimer("bresenhaim::calculatePixelColor");
					pixelColor += (calculatePixelColor(x, y, z, rayDirection, meshIntersectionIndex, innerFaceIntersectionIndex, intersectionPoint) / (float)(supersamplingCoeff*supersamplingCoeff));
					Profiler::finishTimer("bresenhaim::calculatePixelColor");
				}
			}
			pixels[h*imgWidth*4 + w*4] = (unsigned char) (pixelColor.r * 255.0);
			pixels[h*imgWidth*4 + w*4 + 1] = (unsigned char) (pixelColor.g * 255.0);
			pixels[h*imgWidth*4 + w*4 + 2] = (unsigned char) (pixelColor.b * 255.0);
			Profiler::finishTimer("bresenhaim::timePerPixel");
		}
	}
	MImage img;
	img.setPixels(pixels,imgWidth,imgHeight);
	img.writeToFile(outputFilePath);
	img.release();
	delete [] pixels;
}

bool RayTracer::findStartingVoxelIndeces(const MVector& rayDirection, int& x, int& y, int& z)
{
	if(cameraInSceneBB)
	{
		x = initCameraVoxelX;
		y = initCameraVoxelY;
		z = initCameraVoxelZ;
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
		if( sceneBBPlanes[i].rayIntersection(activeCameraData.eye, rayDirection, curTime, curIntersection)
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
	int dimension = voxelParams.voxelsPerDimension;
	while(true)
	{
		if(!(x >= 0 && x < dimension && y >= 0 && y < dimension && z >= 0 && z < dimension))
		{
			break;
		}
		if(pointInVoxelByDirection(point, voxelsData[flatten3dCubeIndex(dimension,x,y,z)], direction))
		{
			return true;
		}
		incrementIndeces(direction, x, y, z);
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
		x = voxelParams.voxelsPerDimension - 1;
	}
	else if (direction == Y_POS)
	{
		y = voxelParams.voxelsPerDimension - 1;

	}
	else if (direction == Z_POS)
	{
		z = voxelParams.voxelsPerDimension - 1;
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

inline bool RayTracer::pointInVoxelByDirection( const MPoint& closestIntersection,VoxelDataT voxel, AxisDirection direction )
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

void RayTracer::incrementIndeces( AxisDirection uDirection, int& x, int& y, int& z )
{
	switch (uDirection)
	{
	case X_POS:
		++x;
		break;
	case X_NEG:
		--x;
		break;
	case Y_POS:
		++y;
		break;
	case Y_NEG:
		--y;
		break;
	case Z_POS:
		++z;
		break;
	case Z_NEG:
		--z;
		break;
	default:
		break;
	}
}

// The function finds the mesh with it intersects given ray, the inner id of the face in the mesh and the intersection point.
// Also it changes the x,y,z indices to match the voxel where the closest intersection happens.
// Returns true if finds
// Return false if it arrives to the scene bounds and doesn't meet any mesh an some point.
bool RayTracer::closestIntersection(const int dimension,const MPoint& raySource,const MVector& rayDirection, int& x, int& y, int& z , int& meshIndex, int& innerFaceId, MPoint& intersection )
{
	Profiler::startTimer("SELF::closestIntersection");
	MPoint nearInt, farInt;
	AxisDirection nearAxisDir, farAxisDir;
	VoxelDataT voxelData;
	int currMeshIndex, currInnerFaceId;
	MPoint currIntersection;
	for(;x >= 0 && x < dimension && y >= 0 && y < dimension && z >= 0 && z < dimension; incrementIndeces(farAxisDir, x, y, z)) {
		voxelData = voxelsData[flatten3dCubeIndex(dimension, x, y, z)];
		if(!voxelData.v->intersectionsWithRay(raySource, rayDirection, nearInt, nearAxisDir, farInt, farAxisDir)) {
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
		Profiler::finishTimer("SELF::closestIntersection");
		return true;
	}
	Profiler::finishTimer("SELF::closestIntersection");
	return false;

	/*pixels[h*imgWidth*4 + w*4] = 0;
	pixels[h*imgWidth*4 + w*4 + 1] = 125;
	pixels[h*imgWidth*4 + w*4 + 2] = 0;*/
}

bool RayTracer::closestIntersectionInVoxel( MPoint raySource, MVector rayDirection, VoxelDataT &voxelData, int &meshIndex, int &innerFaceId, MPoint &intersection )
{
	Profiler::startTimer("SELF::closestIntersectionInVoxel");
	bool res = false;
	double minTime = DBL_MAX;
	double time;
	MPoint triangleVertices[3];
	MIntArray vertexIds;
	int currentMeshIndex, currentFaceIndex;
	MPoint curIntersection;
	for(map<int, vector<int>>::iterator it = voxelData.meshIdToFaceIds.begin(); it != voxelData.meshIdToFaceIds.end(); ++it )
	{
		currentMeshIndex = it->first;
		MFnMesh mesh(meshesData[it->first].dagPath);
		vector<int>& faceIds = it->second;
		for(currentFaceIndex = (int) faceIds.size() - 1; currentFaceIndex >= 0; --currentFaceIndex)
		{
			mesh.getPolygonVertices(faceIds[currentFaceIndex], vertexIds);
			if(vertexIds.length() != 3) {
				continue;
			}
			for (int vi = 0; vi < 3; ++vi) {
				mesh.getPoint(vertexIds[vi], triangleVertices[vi], MSpace::kWorld);
			}
			if(!rayIntersectsTriangle(raySource, rayDirection, triangleVertices, time, curIntersection)
				|| !isPointInVolume(curIntersection, voxelData.v->Min(), voxelData.v->Max()))
			{
				continue;
			}
			if(time < minTime) {
				meshIndex = currentMeshIndex;
				innerFaceId = faceIds[currentFaceIndex];
				intersection = curIntersection;
				minTime = time;
				res = true;
			}
		}
	}
	Profiler::finishTimer("SELF::closestIntersectionInVoxel");
	return res;
}

MColor RayTracer::calculatePixelColor(const int x, const int y, const int z,const MVector& rayDirection, const int meshIndex,const int innerFaceId,const MPoint& intersection)
{
	
	MFnMesh mesh(meshesData[meshIndex].dagPath);
	MIntArray vertexIds;
	MPoint triangleVertices[3];
	MVector triangleNormals[3];
	float us[3];
	float vs[3];

	Profiler::startTimer("calculatePixelColor::gatherMeshGeometry");
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
	}
	Profiler::finishTimer("calculatePixelColor::gatherMeshGeometry");

	double baricentricCoords[3];
	Profiler::startTimer("calculatePixelColor::caclulateBaricentricCoordinates");
	caclulateBaricentricCoordinates(triangleVertices, intersection, baricentricCoords );
	Profiler::finishTimer("calculatePixelColor::caclulateBaricentricCoordinates");

	Profiler::startTimer("calculatePixelColor::gatherMeshMaterial");
	MColor	specularMaterialColor	=	meshesData[meshIndex].specular;
	MColor	ambientMaterialColor	=	meshesData[meshIndex].ambient;
	float	specularPower			=	meshesData[meshIndex].specularPower;
	MColor	diffuseMaterialColor;
	if (meshesData[meshIndex].hasTexture) 
	{
		// get texture color at point using u,v and bilinear filter
		double u = 0, v = 0;
		for (int z = 0; z < 3; ++z)
		{
			u += baricentricCoords[z] * us[z];
			v += baricentricCoords[z] * vs[z];
		}
		diffuseMaterialColor = textureNearesNeighborAtPoint(meshesData[meshIndex].texture, u, v); 
	}
	else 
	{
		diffuseMaterialColor = meshesData[meshIndex].diffuse;
	}
	Profiler::finishTimer("calculatePixelColor::gatherMeshMaterial");

	MVector normalAtPoint = (baricentricCoords[0] * triangleNormals[0] + baricentricCoords[1] * triangleNormals[1] + baricentricCoords[2] * triangleNormals[2]).normal();
	MColor pixelColor = MColor(0,0,0,1);

	LightDataT currLight;
	int currX, currY, currZ;
	MPoint secondIntersection;
	int secondIntersectionMeshIndex, secondIntersectionFaceId;
	for (int li = (int) lightingData.size()- 1; li >= 0; --li)
	{
		MColor currColorComponent(0,0,0,0);
		
		currLight = lightingData[li];
		MColor mixedDiffuse = diffuseMaterialColor * currLight.color * currLight.intencity;
		MColor mixedSpecular = specularMaterialColor * currLight.color * currLight.intencity;	
		currX = x;
		currY = y;
		currZ = z;
		if(currLight.type == LightDataT::AMBIENT)
		{
			MColor mixedAmbient = ambientMaterialColor * currLight.color * currLight.intencity;
			pixelColor = sumColors(mixedAmbient, pixelColor);
		}
		else if(currLight.type == LightDataT::DIRECTIONAL)
		{
			Profiler::startTimer("calculatePixelColor::closestIntersection");
			if(! closestIntersection(voxelParams.voxelsPerDimension, intersection, -currLight.direction, currX, currY, currZ, secondIntersectionMeshIndex, secondIntersectionFaceId, secondIntersection )){
				pixelColor = sumColors(calculateSpecularAndDiffuse(rayDirection, currLight.direction, normalAtPoint, mixedDiffuse, mixedSpecular, specularPower), pixelColor);
			}
			Profiler::finishTimer("calculatePixelColor::closestIntersection");
		}
		else if( currLight.type == LightDataT::POINT)
		{
			MVector lightDirection = intersection - currLight.position;
			MVector lightDirectionNormalized = lightDirection.normal();
			Profiler::startTimer("calculatePixelColor::closestIntersection");
			if(closestIntersection(voxelParams.voxelsPerDimension, intersection, -lightDirectionNormalized, currX, currY, currZ, secondIntersectionMeshIndex, secondIntersectionFaceId, secondIntersection )){
				if((secondIntersection - intersection).length() >= lightDirection.length()) // intersection after the light
				{
					pixelColor = sumColors(calculateSpecularAndDiffuse(rayDirection, lightDirectionNormalized, normalAtPoint, mixedDiffuse, mixedSpecular, specularPower), pixelColor) ;
				}
			}
			else 
			{
				pixelColor = sumColors(calculateSpecularAndDiffuse(rayDirection, lightDirectionNormalized, normalAtPoint, mixedDiffuse, mixedSpecular, specularPower), pixelColor);
			}
			Profiler::finishTimer("calculatePixelColor::closestIntersection");
		}
	}

	return pixelColor;
}

MColor RayTracer::calculateSpecularAndDiffuse(const MVector& viewDirection, MVector lightDirection,  MVector normalAtPoint, MColor mixedDiffuse, MColor mixedSpecular, float specularPower )
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
		MVector reflected = reflectedRay(lightDirection, normalAtPoint);
		k = -(float)(reflected * viewDirection);
		if(k > 0) {
			currColorComponent = sumColors(pow(k, specularPower) * mixedSpecular, currColorComponent);
		}	
	}
	
	return currColorComponent;
}



#pragma region
/*


void printMeshPoints();
void printCamerasInfo();
void printObjectTypesInScene2();
void printObjectTypesInScene();
MPoint getObjectSpaceCentroid(MObject obj);
void getCameraInfo();
void goOverRays();
void calculateSceneBoundingBox();


void RayTracer::triangulateMeshes()
{
MStatus status;
MItDag dagIterator(MItDag::kDepthFirst, MFn::kMesh , &status);
if(badMStatus(status, "MItDag constructor")) { return; }

for(; !dagIterator.isDone(); dagIterator.next())
{
MDagPath dagPath;
status = dagIterator.getPath(dagPath);
if ( badMStatus(status,"MItDag::getPath")) { continue; }

MFnDagNode dagNode(dagPath, &status);
if ( badMStatus(status,"MFnDagNode constructor")) { continue; }

if(dagPath.hasFn(MFn::kMesh))
{
MFnMesh mesh(dagPath, &status);
if ( badMStatus(status,"MFnMesh constructor")) { continue; }
triangulateMesh(mesh);
}
}





}


MMatrix getTransformation(MItDag & it)
{
MStatus mystat;
MDagPath dagPath;
it.getPath(dagPath);
MObject transformNode = dagPath.transform();
MFnDagNode transform(transformNode, &mystat);
MTransformationMatrix matrix(transform.transformationMatrix());
return matrix.asMatrix();
}

void RayTracer::printObjectTypesInScene2()
{
MStatus mystat;

MFnCamera camCreator;
camCreator.setIsOrtho(true);
camCreator.setName("dima_yasha");
MObject theCam = camCreator.create();
MFnCamera camFun(theCam);
MFloatMatrix camProjection = camFun.projectionMatrix();

MPoint ptsToDraw[100];



MMatrix camMatrix;

MItDag camIter(MItDag::kDepthFirst, MFn::kCamera);
while (!camIter.isDone())
{
MFnCamera curCam (camIter.currentItem());
if (curCam.name() != "dima_yasha") {
MGlobal::displayInfo("Not dima_yasha");
camIter.next();
}
MDagPath dagPath;
camIter.getPath(dagPath);
MObject transformNode = dagPath.transform();
MFnDagNode transform(transformNode, &mystat);
MTransformationMatrix matrix(transform.transformationMatrix());
camMatrix = matrix.asMatrix();
break;
}


int i = 0;
MItDag dagIter(MItDag::kDepthFirst, MFn::kMesh);
while (!dagIter.isDone())
{
MDagPath dagPath;
mystat = dagIter.getPath(dagPath);
MObject transformNode = dagPath.transform(&mystat);
if (!mystat && mystat.statusCode() == MStatus::kInvalidParameter) {
MGlobal::displayInfo("no transformNode");
return;
}
MFnDagNode transform(transformNode, &mystat);
if (!mystat) {
MGlobal::displayInfo("MFnDagNode constructor");
return;
}
MTransformationMatrix matrix(transform.transformationMatrix());

MString str;
MMatrix theMatrix = matrix.asMatrix();
MPoint centroid = getObjectSpaceCentroid(dagIter.currentItem());
MPoint inWorldFrame = centroid * theMatrix;
MGlobal::displayInfo(pointToString(inWorldFrame));
ptsToDraw[i++] = (inWorldFrame * camMatrix);
dagIter.next();
}
bool ok = true;
}

void RayTracer::printCamerasInfo()
{

}

MPoint RayTracer::getObjectSpaceCentroid(MObject obj)
{
MStatus mystat;
MItMeshVertex vit(obj, &mystat);

int numVertices = 0;
MPoint aggregatedPoint = MPoint(0,0,0,0);
while (!vit.isDone()) {
MPoint vpoint = vit.position();
aggregatedPoint += vpoint;
numVertices++;
vit.next();
}
if (numVertices != 0) {
aggregatedPoint = aggregatedPoint / ((double)numVertices);
}
aggregatedPoint.w = 1;
return aggregatedPoint;
}


void RayTracer::printObjectTypesInScene()
{
MItDependencyNodes it(MFn::kMesh);
MFnMesh meshFn;
MStatus mystat;
// keep looping until done
MString outstr;


while(!it.isDone())
{
MObject obj = it.item();
MGlobal::displayInfo(obj.apiTypeStr());

MItMeshVertex vit(obj, &mystat);
if (MS::kSuccess != mystat) {
MGlobal::displayInfo("FAIL!");
}

int numVertices = 0;
MPoint aggregatedPoint = MPoint(0,0,0,0);
while (!vit.isDone()) {
MPoint vpoint = vit.position();
outstr.clear();
outstr = "[";
(((((outstr += vpoint.x) += ",") += vpoint.y) += ",") += vpoint.z) += "]";
MGlobal::displayInfo(outstr);
aggregatedPoint += vpoint;
numVertices++;
vit.next();
}

if (numVertices != 0) {
aggregatedPoint = aggregatedPoint / ((double)numVertices);

MFnTransform transFn(obj, &mystat);

if (MS::kSuccess != mystat) {
MGlobal::displayInfo("FAIL!");
}

MTransformationMatrix trans = transFn.transformation();

MMatrix actualTransMat = trans.asMatrix();

outstr.clear();
outstr += "The centroid is in [";
//outstr += res1.x;
//outstr += ",";
//outstr += res1.y;
//outstr += ",";
//outstr += res1.z;
//outstr += ",";
//outstr += res1.w;

outstr += "]";
MGlobal::displayInfo(outstr);
}


// move on to next node
it.next();
}
}

void getLightInfo()
{
MStatus status;
MItDag dagIterator(MItDag::kDepthFirst, MFn::kLight , &status);
if(badMStatus(status, "MItDag constructor"))
{ return; }

vector<MDagPath> lights;

for(; !dagIterator.isDone(); dagIterator.next())
{
MDagPath dagPath;
status = dagIterator.getPath(dagPath);
if ( badMStatus(status,"MItDag::getPath")) { continue; }

MFnDagNode dagNode(dagPath, &status);
if ( badMStatus(status,"MFnDagNode constructor")) { continue; }

if(dagPath.hasFn(MFn::kLight))
{
MFnLight light(dagPath, &status);

if ( badMStatus(status,"MFnLight constructor")) { continue; }
lights.push_back(dagPath);

MColor color = light.color();

if(dagPath.hasFn(MFn::kAmbientLight))
{
cout << "Ambient light " << light.name() << ", color: ["
<< color.r << ", "
<< color.g << ", "
<< color.b << "]\n" << endl;
}
}
}
}

void RayTracer::getCameraInfo()
{
MDagPath cameraPath;
M3dView::active3dView().getCamera( cameraPath );
MStatus status;
MFnCamera camera(cameraPath, &status);
if(badMStatus( status, "MFnCamera c'tor")) { return; }
MVector up = camera.upDirection(MSpace::kWorld);
up.normalize();
MVector view = camera.viewDirection(MSpace::kWorld);
view.normalize();
MVector eye = camera.eyePoint(MSpace::kWorld, &status);
if(badMStatus(status, "MFnCamera.eyePoint")) { return; }
eyePosition = eye;
double focal = camera.focalLength();
double horizontalAperture = camera.horizontalFilmAperture();

rayDirections = new MFloatVector*[imgHeight];
for(int i = 0; i < imgHeight; ++i)
{ rayDirections[i] = new MFloatVector[imgWidth]; }


MVector xAxis = view ^ up;
xAxis.normalize();
MVector yAxis = xAxis ^ view;
yAxis.normalize();
MVector center = eye + view * (focal / 10);
double delta = (horizontalAperture ) / (double)imgWidth;

int halfWidth = imgWidth / 2;
int halfHeight = imgHeight / 2;

MVector leftBottom = center - delta * halfWidth * xAxis - delta * halfHeight * yAxis;
MVector leftBottomCenter = leftBottom + 0.5 * delta * (xAxis + yAxis);

MVector dx = delta * xAxis;
MVector dy = delta * yAxis;


for( int h = 0; h < imgHeight; ++h )
{
for(int w = 0; w < imgWidth; ++w )
{
rayDirections[h][w] = (leftBottomCenter + h * dy + w * dx) - eye;
}
}


MMatrix mat = cameraPath.inclusiveMatrix();
MFloatMatrix cameraMat( mat.matrix );
}

void getPoints(MFnMesh &mesh, int face, int triangle, MPoint *points, MVector * normals )
{
int vertices[3];
mesh.getPolygonTriangleVertices(face, triangle, vertices);
MPointArray allPts;
mesh.getPoints(allPts, MSpace::kWorld);

MColorArray colors;
MColor def(0.5, 0.5, 0.5);

mesh.getColors(colors, NULL, &def);
int len = colors.length();


for(int i = 0; i < 3; ++i)
{
mesh.getPoint(vertices[i], points[i], MSpace::kWorld);
mesh.getVertexNormal(vertices[i],false, normals[i], MSpace::kWorld);

}

}

void RayTracer::goOverRays()
{
unsigned char* pixels = new unsigned char[imgWidth*imgHeight*4];
memset(pixels,0,imgWidth*imgHeight*4);
MStatus status;
MItDag dagIterator(MItDag::kDepthFirst, MFn::kMesh , &status);
if(badMStatus(status, "MItDag constructor")) { return; }

for(; !dagIterator.isDone(); dagIterator.next())
{
MDagPath dagPath;
status = dagIterator.getPath(dagPath);
if ( badMStatus(status,"MItDag::getPath")) { continue; }

MFnDagNode dagNode(dagPath, &status);
if ( badMStatus(status,"MFnDagNode constructor")) { continue; }

if(dagPath.hasFn(MFn::kMesh))
{
MFnMesh mesh(dagPath, &status);
if ( badMStatus(status,"MFnMesh constructor")) { continue; }

MFloatPoint src;
src.x = (float) eyePosition.x;
src.y = (float) eyePosition.y;
src.z = (float) eyePosition.z;



for( int h = 0; h < imgHeight; ++h )
{
for(int w = 0; w < imgWidth; ++w )
{
MFloatPoint intersection;
MFloatVector dir = rayDirections[h][w];
int face;
int triangle;
if (mesh.closestIntersection( src, dir, NULL, NULL, false,MSpace::kWorld, 10000, false, NULL, intersection, NULL, &face, &triangle, NULL, NULL))
{
MPoint pts[3];
MVector normals[3];

getPoints(mesh, face, triangle, pts, normals);

pixels[h*imgWidth*4 + w*4] = 255;
pixels[h*imgWidth*4 + w*4 + 1] = 255;
pixels[h*imgWidth*4 + w*4 + 2] = 255;
}
}
}
}
}

MImage img;
img.setPixels(pixels,imgWidth,imgHeight);
img.writeToFile("C://temp//scene.iff");
img.release();
delete [] pixels;

}

void RayTracer::calculateSceneBoundingBox()
{
MStatus status;
MItDag dagIterator(MItDag::kDepthFirst, MFn::kMesh , &status);
if(badMStatus(status, "MItDag constructor")) { return; }

for(; !dagIterator.isDone(); dagIterator.next())
{
MDagPath dagPath;
status = dagIterator.getPath(dagPath);
if ( badMStatus(status,"MItDag::getPath")) { continue; }

MFnDagNode dagNode(dagPath, &status);
if ( badMStatus(status,"MFnDagNode constructor")) { continue; }

if(dagPath.hasFn(MFn::kMesh))
{
MFnMesh mesh(dagPath, &status);
if ( badMStatus(status,"MFnMesh constructor")) { continue; }

MPointArray pts;
mesh.getPoints(pts, MSpace::kWorld);
for(int pi = pts.length() - 1; pi >= 0; --pi)
{
MPoint current = pts[pi];
if(minScene.x > current.x)
{
minScene.x = current.x;
}
if(minScene.y > current.y)
{
minScene.y = current.y;
}
if(minScene.z > current.z)
{
minScene.z = current.z;
}
if(maxScene.x < current.x)
{
maxScene.x = current.x;
}
if(maxScene.y < current.y)
{
maxScene.y = current.y;
}
if(maxScene.z < current.z)
{
maxScene.z = current.z;
}
}
}
}
}
*/
#pragma endregion

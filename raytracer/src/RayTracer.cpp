#include "RayTracer.h"

#define PRINT_IN_MAYA(arg) MGlobal::displayInfo((arg)); 

#pragma region MyRegion

void* RayTracer::creator() 
{ 
	return new RayTracer; 
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
			/*outstr.clear();
			outstr = "[";
			(((((outstr += vpoint.x) += ",") += vpoint.y) += ",") += vpoint.z) += "]";
			MGlobal::displayInfo(outstr); */
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

#pragma endregion

bool badMStatus(const MStatus& status, const MString& error)
{
	if(status != MStatus::kSuccess)
	{
		status.perror(error);
		return true;
	}
	return false;
}



void triangulateMesh(const MFnMesh& mesh)
{
	MString cmd("polyTriangulate -ch 0 ");
	cmd += mesh.name();
	MGlobal::executeCommand( cmd );
}

void triangulateMeshes()
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

	MSelectionList selected;
	MGlobal::getActiveSelectionList(selected);
	selected.clear();
	MGlobal::setActiveSelectionList(selected);
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
	MVector leftBottomCenter =  leftBottom + 0.5 * delta * (xAxis + yAxis);

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
			src.x = eyePosition.x;
			src.y = eyePosition.y;
			src.z = eyePosition.z;

			

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


MStatus RayTracer::doIt(const MArgList& argList) 
{
	cout << "Running raytracer plugin..." << endl;
	MGlobal::displayInfo("Running raytracer plugin...");


	storeActiveCameraData();
	computeAndStoreImagePlaneData();
	storeLightingData();
	triangulateMeshes();
	computeAndStoreMeshData();
	computeAndStoreSceneBoundingBox();
	voxelizeScene();

	/*
	triangulateMeshes();
	getLightInfo();
	getCameraInfo();
	calculateSceneBoundingBox();
	voxelizeScene();
	goOverRays();
	*/


	MGlobal::displayInfo("Raytracer plugin run finished!");
	cout << "Raytracer plugin run finished!" << endl;  
	return MS::kSuccess;
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

	MPoint centerPoint = activeCameraData.eye + (activeCameraData.focalLengthCm * activeCameraData.viewDir);

	
	double imgAspect = (double)imgWidth / (double)imgHeight;
	//double pixelWidth = (activeCameraData.filmWidthCm)/(double)imgWidth;
	//double pixelHeight = pixelWidth / imgAspect;
	imagePlane.dp = (activeCameraData.filmWidthCm)/(double)imgWidth;

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
	lightingData.push_back(ld);
}

void RayTracer::storeDirectionalLight(MDagPath lightDagPath)
{

	MStatus status;
	MFnLight l(lightDagPath);
	MFloatVector dirFloat = l.lightDirection(lightDagPath.instanceNumber(), MSpace::kWorld, &status);
	CHECK_MSTATUS(status);

	MVector dir;
	dir.x = (double) dirFloat.x;
	dir.y = (double) dirFloat.y;
	dir.z = (double) dirFloat.z;
	//PRINT_IN_MAYA(MString("Directional Light in WF:") + vectorToString(dir));

	LightDataT ld;
	ld.type = LightDataT::DIRECTIONAL;
	ld.color = l.color();
	ld.direction = dir;
	lightingData.push_back(ld);
}

void RayTracer::storePointLight(MDagPath lightDagPath)
{
	MStatus status;
	MFnLight l(lightDagPath);
	MMatrix transform = getDagPathTransformationMatrix(lightDagPath, &status);
	CHECK_MSTATUS(status);
	MPoint pointLightOrigin = MPoint(0,0,0,1);
	MPoint positionWorldFrame = pointLightOrigin * transform;
	//PRINT_IN_MAYA(MString("Point Light in WF:") + pointToString(positionWorldFrame));

	LightDataT ld;
	ld.type = LightDataT::POINT;
	ld.color = l.color();
	ld.position = positionWorldFrame;
	lightingData.push_back(ld);
}

// storeLighting
#pragma endregion 

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
		pair<MPoint,MPoint> boundingBox = computeWfAxisAlignedBoundingBox(dagPath);
		MeshDataT aMesh;
		aMesh.dagPath = dagPath;
		aMesh.min = boundingBox.first;
		aMesh.max = boundingBox.second;
		meshesData.push_back(aMesh);
		PRINT_IN_MAYA(MString("Storing mesh, bb is:") + pointToString(aMesh.min) + "," + pointToString(aMesh.max));
	}
}

void RayTracer::computeAndStoreSceneBoundingBox()
{
	minScene = MPoint( DBL_MAX ,DBL_MAX,DBL_MAX);
	maxScene = MPoint(DBL_MIN, DBL_MIN, DBL_MIN);
	for (int i = 0; i < meshesData.size(); i++) 
	{
		minimize(&(minScene.x), meshesData[i].min.x);
		minimize(&(minScene.y), meshesData[i].min.y);
		minimize(&(minScene.z), meshesData[i].min.z);

		maximize(&(maxScene.x), meshesData[i].max.x);
		maximize(&(maxScene.y), meshesData[i].max.y);
		maximize(&(maxScene.z), meshesData[i].max.z);
	}

	PRINT_IN_MAYA(MString("Scene, bb is:") + pointToString(minScene) + "," + pointToString(maxScene));
}


void RayTracer::voxelizeScene()
{
	computeAndStoreVoxelParams();

}

// TODO: Are voxels supposed to be cubes?
void RayTracer::computeAndStoreVoxelParams()
{
	double sceneSpanX = maxScene.x - minScene.x;
	double sceneSpanY = maxScene.y - minScene.y;
	double sceneSpanZ = maxScene.z - minScene.z;

	maximize(&sceneSpanX, 0.001); // to avoid zero voxel size in case of a flat scene
	maximize(&sceneSpanY, 0.001); // to avoid zero voxel size in case of a flat scene
	maximize(&sceneSpanZ, 0.001); // to avoid zero voxel size in case of a flat scene

	voxelParams.dx = sceneSpanX / voxelParams.voxelsPerDimension;
	voxelParams.dy = sceneSpanY / voxelParams.voxelsPerDimension;
	voxelParams.dz = sceneSpanZ / voxelParams.voxelsPerDimension;

}
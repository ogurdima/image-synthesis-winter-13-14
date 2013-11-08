#include "RayTracer.h"

#include <vector>
using std::vector;

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
		MGlobal::displayInfo(pointToStr(inWorldFrame));
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

MString RayTracer::pointToStr(MPoint p)
{
	MString outstr;
	outstr += "["; 
	outstr += p.x;
	outstr += ",";
	outstr += p.y;
	outstr += ",";
	outstr += p.z;
	outstr += ",";
	outstr += p.w;
	outstr += "]";
	return outstr;
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

	rayDirections = new MFloatVector*[height];
	for(int i = 0; i < height; ++i)
	{ rayDirections[i] = new MFloatVector[width]; }


	MVector xAxis = view ^ up;
	xAxis.normalize();
	MVector yAxis = xAxis ^ view;
	yAxis.normalize();
	MVector center = eye + view * (focal / 10);
	double delta = (horizontalAperture ) / (double)width;

	int halfWidth = width / 2;
	int halfHeight = height / 2;

	MVector leftBottom = center - delta * halfWidth * xAxis - delta * halfHeight * yAxis;
	MVector leftBottomCenter =  leftBottom + 0.5 * delta * (xAxis + yAxis);

	MVector dx = delta * xAxis;
	MVector dy = delta * yAxis;


	for( int h = 0; h < height; ++h )
	{
		for(int w = 0; w < width; ++w )
		{
			rayDirections[h][w] = (leftBottomCenter + h * dy + w * dx) - eye;
		}
	}
	

	MMatrix mat = cameraPath.inclusiveMatrix();
    MFloatMatrix cameraMat( mat.matrix );
}

void RayTracer::goOverRays()
{
	unsigned char* pixels = new unsigned char[width*height*4];
	memset(pixels,0,width*height*4);
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

			for( int h = 0; h < height; ++h )
			{
				for(int w = 0; w < width; ++w )
				{
					MFloatPoint intersection;
					MFloatVector dir = rayDirections[h][w];
					if (mesh.closestIntersection( src, dir, NULL, NULL, false,MSpace::kWorld, 10000, false, NULL, intersection, NULL, NULL, NULL, NULL, NULL))
					{
						pixels[h*width*4 + w*4] = 255;
						pixels[h*width*4 + w*4 + 1] = 255;
						pixels[h*width*4 + w*4 + 2] = 255;
					}
				}
			}
		}
	}

	MImage img;
	img.setPixels(pixels,width,height);
	img.writeToFile("C://temp//scene.iff");
	img.release();
	delete [] pixels;

}

MStatus RayTracer::doIt(const MArgList& argList) 
{
	//MGlobal::displayInfo("Welcome to Image Synthesis course");
	//printCamerasInfo();
	triangulateMeshes();
	getLightInfo();
	getCameraInfo();

	goOverRays();

	MGlobal::displayInfo("DONE!");
	return MS::kSuccess;

}
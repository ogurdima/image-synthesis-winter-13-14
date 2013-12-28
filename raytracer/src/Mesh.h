#pragma once
#include <maya/MPoint.h>
#include <maya/MColor.h>
#include <maya/MPointArray.h>
#include <maya/MVectorArray.h>
#include <maya/MFloatArray.h>
#include <maya/MImage.h>
#include <maya/MDagPath.h>
#include <vector>
#include "Material.h"


using std::vector;


struct Face
	{
		MPointArray vertices;
		MVectorArray normals;
		MFloatArray us;
		MFloatArray vs;
	};


struct MeshDataT
	{
		//MDagPath	dagPath;

		MPoint		max;		// WS axis aligned bounding box min
		MPoint		min;		// WS axis aligned bounding box max

		/*bool		hasTexture;
		MImage*		texture;
		MColor		diffuse;
		MColor		specular;
		MColor		ambient;
		float		specularPower;

		bool		useHalfVector;
		float		eccentricity;*/
		
		Material	material;

		vector<Face> faces;


	};


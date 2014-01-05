#pragma once

#include <maya\MString.h>
#include <maya\MFnDependencyNode.h>
#include <maya\MPlug.h>
#include <maya\MPlugArray.h>
#include <maya\MColor.h>
#include <vector>
#include <maya/MAngle.h>
#include <maya/MFnTransform.h>
#include <maya/MItDag.h>
#include <maya/MFnCamera.h>
#include <maya/MGlobal.h>
#include <maya/MString.h>
#include <maya/MStringArray.h>
#include <maya/MObject.h>
#include <maya/MObjectArray.h>
#include <maya/MIntArray.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MArgList.h>
#include <maya/MStatus.h>
#include <maya/MDagPath.h>
#include <maya/MFnMesh.h>
#include <maya/MFnLambertShader.h>
#include <maya/MFnPhongShader.h>
#include <maya/MFnBlinnShader.h>
#include <maya/MFnIkJoint.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MMatrix.h>
#include <maya/MFnSkinCluster.h>
#include <maya/MItDependencyNodes.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MFloatArray.h>
#include <maya/MDagPathArray.h>
#include <maya/MPointArray.h>
#include <maya/MItGeometry.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MSelectionList.h>
#include <maya/MItSelectionList.h>
#include <maya/MItDependencyGraph.h>
#include <maya/MTime.h>
#include <maya/MImage.h>
#include <maya/MAnimControl.h>
#include <maya/MAnimUtil.h>
#include <maya/MRenderUtil.h>
#include <maya/MQuaternion.h>
#include <maya/MFnBlendShapeDeformer.h>
#include <maya/MBoundingBox.h>
#include <maya/MDagModifier.h>

#define PRECISION 0.0001

//#include "mayaExportLayer.h"



typedef enum {MT_LAMBERT,MT_PHONG} MaterialType;

typedef enum {TOT_REPLACE,TOT_MODULATE,TOT_ADD,TOT_ALPHABLEND} TexOpType;

typedef enum {TAM_CLAMP,TAM_BORDER,TAM_WRAP,TAM_MIRROR} TexAddressMode;

class Texture
{
public:
	//constructor
	Texture() {
		scale_u = scale_v = 1;
		scroll_u = scroll_v = 0;
		rot = 0;
		am_u = am_v = TAM_CLAMP;
	}
	//destructor
	~Texture(){};

	//public members
	MString filename;
	MString absFilename;
	TexOpType opType;
	MString uvsetName;
	int uvsetIndex;
	TexAddressMode am_u,am_v;
	double scale_u,scale_v;
	double scroll_u,scroll_v;
	double rot;
};


/***** Class Material *****/
class Material
{
public:
	//constructor
	Material();
	//destructor
	~Material();
	//get material name
	MString& name();
	//clear material data
	void clear();
	//load material data
	MStatus load(MFnDependencyNode* pShader,MStringArray& uvsets);
	//load a specific material type
	//MStatus loadSurfaceShader(MFnDependencyNode* pShader);
	MStatus loadLambert(MFnDependencyNode* pShader);
	MStatus loadPhong(MFnDependencyNode* pShader);
	//MStatus loadBlinn(MFnDependencyNode* pShader);
	//MStatus loadCgFxShader(MFnDependencyNode* pShader);

	void toDefault();

public:
	//load texture data
	MStatus loadTexture(MFnDependencyNode* pTexNode,TexOpType& opType,MStringArray& uvsets);

	//MString m_name;
	MaterialType type;
	MColor ambient, diffuse, specular, emissive;
	float diffuseCoeff;
	bool isTransparent;
	bool isReflective;
	bool isTextured;

	float cosPower;
	float transparency;
	float refractiveIndex;

	double kr0;

	float reflectivity;

	//bool m_isMultiTextured;
	//std::vector<Texture> m_textures;
	MImage * texture;
};

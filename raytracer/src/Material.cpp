#include "Material.h"



// Constructor
Material::Material()
{
	texture = NULL;
	clear();
}


// Destructor
Material::~Material()
{
}



// Clear data
void Material::clear()
{
	type = MT_LAMBERT;
	isTransparent = false;
	isTextured = false;
	//isMultiTextured = false;
	ambient = MColor(0,0,0,0);
	diffuse = MColor(0,0,0,0);
	specular = MColor(0,0,0,0);
	emissive = MColor(0,0,0,0);
	if(texture){
		delete texture;
		texture = NULL;
	}
}


bool getLambertShaderTexture(MFnDependencyNode* lambert, MImage& img)
{
	MPlugArray plugs;
	lambert->findPlug("color").connectedTo(plugs, true, false);
	for(uint i = 0; i < plugs.length(); i++)
	{
		if (plugs[i].node().hasFn(MFn::kFileTexture))
		{
			img.release();
			img.readFromTextureNode( plugs[i].node() );
			return true;
		}
	}
	return false;
}

// Load material data
MStatus Material::load(MFnDependencyNode* pShader,MStringArray& uvsets)
{
	MStatus stat;
	clear();
	////read material name, adding the requested prefix
	//MString tmpStr = params.matPrefix;
	//if (tmpStr != "")
	//	tmpStr += "/";
	//tmpStr += pShader->name();
	//MStringArray tmpStrArray;
	//tmpStr.split(':',tmpStrArray);
	//name = "";
	//for (int i=0; i<tmpStrArray.length(); i++)
	//{
	//	name += tmpStrArray[i];
	//	if (i < tmpStrArray.length()-1)
	//		name += "_";
	//}

	//check if we want to export with lighting off option
	//lightingOff = params.lightingOff;

	// GET MATERIAL DATA

	// Check material type
	if (pShader->object().hasFn(MFn::kPhong))
	{
		stat = loadPhong(pShader);
	}/*
	 else if (pShader->object().hasFn(MFn::kBlinn))
	 {
	 stat = loadBlinn(pShader);
	 }*/
	else if (pShader->object().hasFn(MFn::kLambert))
	{
		stat = loadLambert(pShader);
	}
	/*else if (pShader->object().hasFn(MFn::kPluginHwShaderNode))
	{
	stat = loadCgFxShader(pShader);
	}
	else
	{
	stat = loadSurfaceShader(pShader);
	}*/



	// Get textures data
#pragma region Old texture reading
	//MPlugArray colorSrcPlugs;
	//MPlugArray texSrcPlugs;
	//MPlugArray placetexSrcPlugs;
	//if (isTextured)
	//{
	//	// Translate multiple textures if material is multitextured
	//	if (isMultiTextured)
	//	{
	//		// Get layered texture node
	//		MFnDependencyNode* pLayeredTexNode = NULL;
	//		if (type == MT_SURFACE_SHADER)
	//			pShader->findPlug("outColor").connectedTo(colorSrcPlugs,true,false);
	//		else
	//			pShader->findPlug("color").connectedTo(colorSrcPlugs,true,false);
	//		for (int i=0; i<colorSrcPlugs.length(); i++)
	//		{
	//			if (colorSrcPlugs[i].node().hasFn(MFn::kLayeredTexture))
	//			{
	//				pLayeredTexNode = new MFnDependencyNode(colorSrcPlugs[i].node());
	//				continue;
	//			}
	//		}

	//		// Get inputs to layered texture
	//		MPlug inputsPlug = pLayeredTexNode->findPlug("inputs");

	//		// Scan inputs and export textures
	//		for (int i=inputsPlug.numElements()-1; i>=0; i--)
	//		{
	//			MFnDependencyNode* pTextureNode = NULL;
	//			// Search for a connected texture
	//			inputsPlug[i].child(0).connectedTo(colorSrcPlugs,true,false);
	//			for (int j=0; j<colorSrcPlugs.length(); j++)
	//			{
	//				if (colorSrcPlugs[j].node().hasFn(MFn::kFileTexture))
	//				{
	//					pTextureNode = new MFnDependencyNode(colorSrcPlugs[j].node());
	//					continue;
	//				}
	//			}

	//			// Translate the texture if it was found
	//			if (pTextureNode)
	//			{
	//				// Get blend mode
	//				TexOpType opType;
	//				short bm;
	//				inputsPlug[i].child(2).getValue(bm);
	//				switch(bm)
	//				{				
	//				case 0:
	//					opType = TOT_REPLACE;
	//					break;
	//				case 1:
	//					opType = TOT_ALPHABLEND;
	//					break;				
	//				case 4:
	//					opType = TOT_ADD;
	//					break;
	//				case 6:
	//					opType = TOT_MODULATE;
	//					break;
	//				default:
	//					opType = TOT_MODULATE;
	//				}

	//				stat = loadTexture(pTextureNode,opType,uvsets);
	//				delete pTextureNode;
	//				if (MS::kSuccess != stat)
	//				{
	//					std::cout << "Error loading layered texture\n";
	//					std::cout.flush();
	//					delete pLayeredTexNode;
	//					return MS::kFailure;
	//				}
	//			}
	//		}
	//		if (pLayeredTexNode)
	//			delete pLayeredTexNode;
	//	}
	//	// Else translate the single texture
	//	else
	//	{
	//		// Get texture node
	//		MFnDependencyNode* pTextureNode = NULL;
	//		if (type == MT_SURFACE_SHADER)
	//			pShader->findPlug("outColor").connectedTo(colorSrcPlugs,true,false);
	//		else
	//			pShader->findPlug("color").connectedTo(colorSrcPlugs,true,false);
	//		for (int i=0; i<colorSrcPlugs.length(); i++)
	//		{
	//			if (colorSrcPlugs[i].node().hasFn(MFn::kFileTexture))
	//			{
	//				pTextureNode = new MFnDependencyNode(colorSrcPlugs[i].node());
	//				continue;
	//			}
	//		}
	//		if (pTextureNode)
	//		{
	//			TexOpType opType = TOT_MODULATE;
	//			stat = loadTexture(pTextureNode,opType,uvsets);
	//			delete pTextureNode;
	//			if (MS::kSuccess != stat)
	//			{
	//				std::cout << "Error loading texture\n";
	//				std::cout.flush();
	//				return MS::kFailure;
	//			}
	//		}
	//	}  
	// }
#pragma endregion



	return MS::kSuccess;
}


// Load a surface shader
//MStatus Material::loadSurfaceShader(MFnDependencyNode *pShader)
//{
//	type = MT_SURFACE_SHADER;
//	MPlugArray colorSrcPlugs;
//	// Check if material is textured
//	pShader->findPlug("outColor").connectedTo(colorSrcPlugs,true,false);
//	for (int i=0; i<colorSrcPlugs.length(); i++)
//	{
//		if (colorSrcPlugs[i].node().hasFn(MFn::kFileTexture))
//		{
//			isTextured = true;
//			continue;
//		}
//		else if (colorSrcPlugs[i].node().hasFn(MFn::kLayeredTexture))
//		{
//			isTextured = true;
//			isMultiTextured = true;
//			continue;
//		}
//	}

//	// Check if material is transparent
//	float trasp;
//	pShader->findPlug("outTransparencyR").getValue(trasp);
//	if (pShader->findPlug("outTransparency").isConnected() || trasp>0.0f)
//		isTransparent = true;

//	// Get material colours
//	if (isTextured)
//		diffuse = MColor(1.0,1.0,1.0,1.0);
//	else
//	{
//		pShader->findPlug("outColorR").getValue(diffuse.r);
//		pShader->findPlug("outColorG").getValue(diffuse.g);
//		pShader->findPlug("outColorB").getValue(diffuse.b);
//		float trasp;
//		pShader->findPlug("outTransparencyR").getValue(trasp);
//		diffuse.a = 1.0 - trasp;
//	}
//	ambient = MColor(0,0,0,1);
//	emissive = MColor(0,0,0,1);
//	specular = MColor(0,0,0,1);
//	return MS::kSuccess;
//}

// Load a lambert shader
MStatus Material::loadLambert(MFnDependencyNode *pShader)
{
	MPlugArray colorSrcPlugs;
	type = MT_LAMBERT;
	MFnLambertShader* pLambert = new MFnLambertShader(pShader->object());
	// Check if material is textured

	MImage* textureImage = new MImage();
	if (getLambertShaderTexture(pShader, *textureImage)) 
	{
		// store the texture image
		texture = textureImage;
		isTextured = true;
	}
	else 
	{
		// assign diffuse and ambient color
		isTextured = false;
		delete textureImage;
		diffuse = pLambert->color();
	}

	diffuseCoeff = pLambert->diffuseCoeff();
	//ambient colour
	ambient = pLambert->ambientColor();
	//emissive colour
	emissive = pLambert->incandescence();
	//specular colour
	specular = MColor(0,0,0,1);


	// Check if material is transparent
	if (pLambert->findPlug("transparency").isConnected() || pLambert->transparency().r>0.0f){
		isTransparent = true;
		transparancy = pLambert->transparency().r;
	}

	delete pLambert;
	return MS::kSuccess;
}

void Material::toDefault()
{
	clear();
	type = MT_LAMBERT;
	ambient = MColor(0.1f, 0.1f, 0.1f);
	diffuse = MColor(0.f, 0.f, 0.8f);
	specular = MColor(0.9f,0.9f,0.9f);
	emissive = MColor(0.f,0.f,0.f);
	cosPower = 10;
	isTextured = false;
	isTransparent = false;
}

	// Load a phong shader
MStatus Material::loadPhong(MFnDependencyNode *pShader)
{
	MPlugArray colorSrcPlugs;
	type = MT_PHONG;
	MFnPhongShader* pPhong = new MFnPhongShader(pShader->object());
	// Check if material is textured
	MImage* textureImage = new MImage();
	if (getLambertShaderTexture(pShader, *textureImage)) 
	{
		// store the texture image
		texture = textureImage;
		isTextured = true;
	}
	else 
	{
		// assign diffuse and ambient color
		isTextured = false;
		delete textureImage;
		diffuse = pPhong->color();
	}

	diffuseCoeff = pPhong->diffuseCoeff();

	//ambient colour
	ambient = pPhong->ambientColor();
	//emissive colour
	emissive = pPhong->incandescence();
	//specular colour
	specular = pPhong->specularColor();

	cosPower = pPhong->cosPower();

	// Check if material is transparent
	if (pPhong->findPlug("transparency").isConnected() || pPhong->transparency().r>0.0f){
		isTransparent = true;
		transparancy = pPhong->transparency().r;
	}

	delete pPhong;
	return MS::kSuccess;
}

// load a blinn shader
//MStatus Material::loadBlinn(MFnDependencyNode *pShader)
//{
//	MPlugArray colorSrcPlugs;
//	type = MT_BLINN;
//	MFnBlinnShader* pBlinn = new MFnBlinnShader(pShader->object());
//	// Check if material is textured
//	pBlinn->findPlug("color").connectedTo(colorSrcPlugs,true,false);
//	for (int i=0; i<colorSrcPlugs.length(); i++)
//	{
//		if (colorSrcPlugs[i].node().hasFn(MFn::kFileTexture))
//		{
//			isTextured = true;
//			continue;
//		}
//		/*else if (colorSrcPlugs[i].node().hasFn(MFn::kLayeredTexture))
//		{
//			isTextured = true;
//			isMultiTextured = true;
//			continue;
//		}*/
//	}

//	// Check if material is transparent
//	if (pBlinn->findPlug("transparency").isConnected() || pBlinn->transparency().r>0.0f)
//		isTransparent = true;

//	// Get material colours
//	//diffuse colour
//	if (isTextured)
//		diffuse = MColor(1.0,1.0,1.0,1.0);
//	else
//	{
//		diffuse = pBlinn->color();
//		diffuse.a = 1.0 - pBlinn->transparency().r;
//	}
//	//ambient colour
//	ambient = pBlinn->ambientColor();
//	//emissive colour
//	emissive = pBlinn->incandescence();
//	//specular colour
//	specular = pBlinn->specularColor();
//	specular.a = 128.0-(128.0*pBlinn->eccentricity());
//	delete pBlinn;
//	return MS::kSuccess;
//}

// load a cgFx shader
//MStatus Material::loadCgFxShader(MFnDependencyNode *pShader)
//{
//	type = MT_CGFX;
//	// Create a default white lambert
//	isTextured = false;
//	isMultiTextured = false;
//	diffuse = MColor(1.0,1.0,1.0,1.0);
//	specular = MColor(0,0,0,1);
//	emissive = MColor(0,0,0,1);
//	ambient = MColor(0,0,0,1);
//	return MS::kSuccess;
//}

// Load texture data from a texture node
//MStatus Material::loadTexture(MFnDependencyNode* pTexNode,TexOpType& opType,MStringArray& uvsets)
//{
//	Texture tex;
//	// Get texture filename
//	MString filename, absFilename;
//	MRenderUtil::exactFileTextureName(pTexNode->object(),absFilename);
//	filename = absFilename.substring(absFilename.rindex('/')+1,absFilename.length()-1);
//	MString command = "toNativePath(\"";
//	command += absFilename;
//	command += "\")";
//	MGlobal::executeCommand(command,absFilename);
//	tex.absFilename = absFilename;
//	tex.filename = filename;
//	tex.uvsetIndex = 0;
//	tex.uvsetName = "";
//	// Set texture operation type
//	tex.opType = opType;
//	// Get connections to uvCoord attribute of texture node
//	MPlugArray texSrcPlugs;
//	pTexNode->findPlug("uvCoord").connectedTo(texSrcPlugs,true,false);
//	// Get place2dtexture node (if connected)
//	MFnDependencyNode* pPlace2dTexNode = NULL;
//	for (int j=0; j<texSrcPlugs.length(); j++)
//	{
//		if (texSrcPlugs[j].node().hasFn(MFn::kPlace2dTexture))
//		{
//			pPlace2dTexNode = new MFnDependencyNode(texSrcPlugs[j].node());
//			continue;
//		}
//	}
//	// Get uvChooser node (if connected)
//	MFnDependencyNode* pUvChooserNode = NULL;
//	if (pPlace2dTexNode)
//	{
//		MPlugArray placetexSrcPlugs;
//		pPlace2dTexNode->findPlug("uvCoord").connectedTo(placetexSrcPlugs,true,false);
//		for (int j=0; j<placetexSrcPlugs.length(); j++)
//		{
//			if (placetexSrcPlugs[j].node().hasFn(MFn::kUvChooser))
//			{
//				pUvChooserNode = new MFnDependencyNode(placetexSrcPlugs[j].node());
//				continue;
//			}
//		}
//	}
//	// Get uvset index
//	if (pUvChooserNode)
//	{
//		bool foundMesh = false;
//		bool foundUvset = false;
//		MPlug uvsetsPlug = pUvChooserNode->findPlug("uvSets");
//		MPlugArray uvsetsSrcPlugs;
//		for (int i=0; i<uvsetsPlug.evaluateNumElements() && !foundMesh; i++)
//		{
//			uvsetsPlug[i].connectedTo(uvsetsSrcPlugs,true,false);
//			for (int j=0; j<uvsetsSrcPlugs.length() && !foundMesh; j++)
//			{
//				if (uvsetsSrcPlugs[j].node().hasFn(MFn::kMesh))
//				{
//					uvsetsSrcPlugs[j].getValue(tex.uvsetName);
//					for (int k=0; k<uvsets.length() && !foundUvset; k++)
//					{
//						if (uvsets[k] == tex.uvsetName)
//						{
//							tex.uvsetIndex = k;
//							foundUvset = true;
//						}
//					}
//				}
//			}
//		}
//	}
//	// Get texture options from Place2dTexture node
//	if (pPlace2dTexNode)
//	{
//		// Get address mode
//		//U
//		bool wrapU, mirrorU;
//		pPlace2dTexNode->findPlug("wrapU").getValue(wrapU);
//		pPlace2dTexNode->findPlug("mirrorU").getValue(mirrorU);
//		if (mirrorU)
//			tex.au = TAMIRROR;
//		else if (wrapU)
//			tex.au = TAWRAP;
//		else
//			tex.au = TACLAMP;
//		// V
//		bool wrapV,mirrorV;
//		pPlace2dTexNode->findPlug("wrapV").getValue(wrapV);
//		pPlace2dTexNode->findPlug("mirrorV").getValue(mirrorV);
//		if (mirrorV)
//			tex.av = TAMIRROR;
//		else if (wrapV)
//			tex.av = TAWRAP;
//		else
//			tex.av = TACLAMP;
//		// Get texture scale
//		double covU,covV;
//		pPlace2dTexNode->findPlug("coverageU").getValue(covU);
//		pPlace2dTexNode->findPlug("coverageV").getValue(covV);
//		tex.scale_u = covU;
//		if (fabs(tex.scale_u) < PRECISION)
//			tex.scale_u = 0;
//		tex.scale_v = covV;
//		if (fabs(tex.scale_v) < PRECISION)
//			tex.scale_v = 0;
//		// Get texture scroll
//		double transU,transV;
//		pPlace2dTexNode->findPlug("translateFrameU").getValue(transU);
//		pPlace2dTexNode->findPlug("translateFrameV").getValue(transV);
//		tex.scroll_u = -0.5 * (covU-1.0)/covU - transU/covU;
//		if (fabs(tex.scroll_u) < PRECISION)
//			tex.scroll_u = 0;
//		tex.scroll_v = 0.5 * (covV-1.0)/covV + transV/covV;
//		if (fabs(tex.scroll_v) < PRECISION)
//			tex.scroll_v = 0;
//		// Get texture rotation
//		double rot;
//		pPlace2dTexNode->findPlug("rotateFrame").getValue(rot);
//		tex.rot = -rot;
//		if (fabs(rot) < PRECISION)
//			tex.rot = 0;
//	}
//	// add texture to material texture list
//	textures.push_back(tex);
//	// free up memory
//	if (pUvChooserNode)
//		delete pUvChooserNode;
//	if (pPlace2dTexNode)
//		delete pPlace2dTexNode;

//	return MS::kSuccess;
//}


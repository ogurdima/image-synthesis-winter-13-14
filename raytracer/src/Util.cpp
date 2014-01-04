#include "Util.h"

#include <math.h>
#include <stdio.h>

namespace util 
{

	MString pointToString(MPoint p)
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

	MString vectorToString(MVector p)
	{
		MString outstr;
		outstr += "["; 
		outstr += p.x;
		outstr += ",";
		outstr += p.y;
		outstr += ",";
		outstr += p.z;
		outstr += "]";
		return outstr;
	}

	MString colorToString(MColor c)
	{
		MString outstr;
		outstr += "["; 
		outstr += c.r;
		outstr += ",";
		outstr += c.g;
		outstr += ",";
		outstr += c.b;
		outstr += ",";
		outstr += c.a;
		outstr += "]";
		return outstr;
	}

	MMatrix getDagPathTransformationMatrix(MDagPath dagPath, MStatus* statusPtr)
	{
		MStatus status;
		MObject transformNode = dagPath.transform(&status);
		if (status != MS::kSuccess) {
			if (statusPtr != NULL) {
				*statusPtr = status;
			}
			return MMatrix();
		}
		MFnDagNode transform(transformNode, &status);
		if (status != MS::kSuccess) {
			if (statusPtr != NULL) {
				*statusPtr = status;
			}
			return MMatrix();
		}
		MTransformationMatrix matrix(transform.transformationMatrix(&status));
		if (status != MS::kSuccess) {
			if (statusPtr != NULL) {
				*statusPtr = status;
			}
			return MMatrix();
		}
		return matrix.asMatrix();
	}

	pair<MPoint, MPoint> computeWfAxisAlignedBoundingBox(MDagPath meshPath, MStatus* statusPtr)
	{
		MStatus status;
		MFnMesh mesh(meshPath, &status);
		if ( status != MS::kSuccess) {
			if (statusPtr != NULL) {
				*statusPtr = status;
			}
			return pair<MPoint,MPoint>(MPoint(), MPoint());
		}

		MPointArray pts;
		MPoint min = MPoint( DBL_MAX ,DBL_MAX, DBL_MAX);
		MPoint max = MPoint( -DBL_MAX ,-DBL_MAX, -DBL_MAX);
		mesh.getPoints(pts, MSpace::kWorld);
		for(int pi = pts.length() - 1; pi >= 0; --pi)
		{
			MPoint current = pts[pi];
			minimize(&(min.x), current.x);
			minimize(&(min.y), current.y);
			minimize(&(min.z), current.z);

			maximize(&(max.x), current.x);
			maximize(&(max.y), current.y);
			maximize(&(max.z), current.z);
		}
		return pair<MPoint,MPoint>(min, max);
	}

	double minimize(double* oldPtr, double newVal)
	{
		if (oldPtr == NULL) {
			return DBL_MAX;
		}
		if (*oldPtr > newVal) {
			*oldPtr = newVal;
		}
		return *oldPtr;
	}

	double maximize(double* oldPtr, double newVal)
	{
		if (oldPtr == NULL) {
			return -DBL_MAX;
		}
		if (*oldPtr < newVal) {
			*oldPtr = newVal;
		}
		return *oldPtr;
	}

	bool intervalsOverlap(double x1, double y1, double x2, double y2)
	{
		double l1, r1, l2, r2;
		if (x1 < y1) {
			l1 = x1;
			r1 = y1;
		}
		else {
			l1 = y1;
			r1 = x1;
		}

		if (x2 < y2) {
			l2 = x2;
			r2 = y2;
		}
		else {
			l2 = y2;
			r2 = x2;
		}


		if (r1 < (l2 - DOUBLE_NUMERICAL_THRESHHOLD) || r2 < (l1 - DOUBLE_NUMERICAL_THRESHHOLD)) {
			return false;
		}
		return true;			
	}

	bool isPointInVolume(const MPoint& point, const MPoint& minVolume, const MPoint& maxVolume)
	{
		for(int i = 0; i < 3; ++i)
		{
			if( point[i] < minVolume[i] - DOUBLE_NUMERICAL_THRESHHOLD || point[i] > maxVolume[i] + DOUBLE_NUMERICAL_THRESHHOLD)
				return false;
		}
		return true;
	}

	bool valueInInterval( double value, double intervalMin, double intervalMax )
	{
		return !( value < intervalMin || value > intervalMax);
	}

	bool pointInRectangle( AxisDirection projectionDirection, const MPoint& point, const MPoint& minPoint, const MPoint& maxPoint )
	{
		switch (projectionDirection)
		{
		case X_NEG:
		case X_POS:
			return valueInInterval(point.y, minPoint.y, maxPoint.y) && valueInInterval(point.z, minPoint.z, maxPoint.z);
		case Y_NEG:
		case Y_POS:
			return valueInInterval(point.x, minPoint.x, maxPoint.x) && valueInInterval(point.z, minPoint.z, maxPoint.z);
		case Z_NEG:
		case Z_POS:
			return valueInInterval(point.y, minPoint.y, maxPoint.y) && valueInInterval(point.x, minPoint.x, maxPoint.x);
		default:
			break;
		}
		return false;
	}

	MColor textureNearesNeighborAtPoint(const MImage* texture, double u, double v, bool repeat)
	{
		MColor res = MColor(0,0,0);
		unsigned char * pixs = texture->pixels();
		unsigned int width, height;
		texture->getSize(width, height);

		double fixedU, fixedV;
		if (repeat) {
			fixedU = fmod(u, 1.0);
			fixedV = fmod(v, 1.0);
		}
		else {
			fixedU = (u < 0) ? 0 : (u > 1) ? 1 : u;
			fixedV = (v < 0) ? 0 : (v > 1) ? 1 : v;
		}

		if(texture->pixelType() == MImage::MPixelType::kByte )
		{
			int x = (int) ( ((double)(width - 1) * fixedU) + 0.49 );
			int y = (int) ( ((double)(height - 1) * fixedV) + 0.49 );
			
			res.r = (float) ((float)pixs[(y * width + x)*4] / 255.0);
			res.g = (float) ((float)pixs[(y * width + x)*4 + 1] / 255.0);
			res.b = (float) ((float)pixs[(y * width + x)*4 + 2]/ 255.0);
		}
		return res;
	}

	MColor getBilinearFilteredPixelColor( const MImage* texture, double u, double v )
	{
		unsigned int w, h;
		MColor res = MColor(0,0,0);
		unsigned char * pixs = texture->pixels();
		texture->getSize(w, h);
		u = u * (w - 1) - 0.5;
		v = v * (h - 1) - 0.5;
		int x = max((int) floor(u), 0);
		int y = max((int) floor(v), 0);
		int xNext = min(x + 1, (int)w - 1);
		int yNext = min(y + 1, (int)h - 1);
		double u_ratio = u - x;
		double v_ratio = v - y; 
		double u_opposite = 1 - u_ratio;
		double v_opposite = 1 - v_ratio;
		for (int i = 0; i < 3; ++i)
		{
			res[i] = (float)(( ((float)pixs[(y * w + x)*4 +i] / 255.0)   * u_opposite  + ((float)pixs[(y * w + xNext)*4 + i] / 255.0)   * u_ratio) * v_opposite +
				(((float)pixs[((yNext) * w + x)*4 + i] / 255.0) * u_opposite  + ((float)pixs[(yNext * w + xNext)*4 + i] / 255.0) * u_ratio) * v_ratio);
		}
		
		return res;
	}

	bool rayIntersectsTriangle(const MPoint& raySrc,const MVector& rayDirection, const MPoint triangleVertices[3], double& time, MPoint& intersection) 
	{
		double a,f,u,v;
		MVector edge01(triangleVertices[1] - triangleVertices[0]);
		MVector edge02(triangleVertices[2] - triangleVertices[0]);
		MVector h = rayDirection ^ edge02;
		a = edge01 * h;
		if (abs(a) < DOUBLE_NUMERICAL_THRESHHOLD) {
			return(false);
		}
		f = 1/a;
		MVector s = raySrc - triangleVertices[0];
		u = f * (s * h);

		if (u < 0.0 || u > 1.0) {
			return(false);
		}
		MVector q = s ^ edge01;
		v = f * (rayDirection * q);
		if (v < 0.0 || u + v > 1.0) {
			return(false);
		}
		// at this stage we can compute t to find out where
		// the intersection point is on the line
		//double t = f * innerProduct(e2,q);
		double localTime = f * (edge02 * q);
		if (localTime > DOUBLE_NUMERICAL_THRESHHOLD) // ray intersection
		{
			time = localTime;
			intersection = raySrc + time * rayDirection;
			return true;
		}
		else {
			// this means that there is a line intersection
			// but not a ray intersection
			return false;
		}
	}

	bool rayIntersectsTriangle(const MPoint& raySrc,const MVector& rayDirection, const MPointArray& triangleVertices, double& time, MPoint& intersection) 
	{
		double a,f,u,v;
		MVector edge01(triangleVertices[1] - triangleVertices[0]);
		MVector edge02(triangleVertices[2] - triangleVertices[0]);
		MVector h = rayDirection ^ edge02;
		a = edge01 * h;
		if (abs(a) < DOUBLE_NUMERICAL_THRESHHOLD) {
			return(false);
		}
		f = 1/a;
		MVector s = raySrc - triangleVertices[0];
		u = f * (s * h);

		if (u < 0.0 || u > 1.0) {
			return(false);
		}
		MVector q = s ^ edge01;
		v = f * (rayDirection * q);
		if (v < 0.0 || u + v > 1.0) {
			return(false);
		}
		// at this stage we can compute t to find out where
		// the intersection point is on the line
		//double t = f * innerProduct(e2,q);
		double localTime = f * (edge02 * q);
		if (localTime > DOUBLE_NUMERICAL_THRESHHOLD) // ray intersection
		{
			time = localTime;
			intersection = raySrc + time * rayDirection;
			return true;
		}
		else {
			// this means that there is a line intersection
			// but not a ray intersection
			return false;
		}
	}
		
	void calculateBaricentricCoordinates(const MPointArray& triangleVertices, const MPoint& point, double baricentricCoords[3] )
	{
		
		MVector e01 = (triangleVertices[1] - triangleVertices[0]);
		MVector e02 = (triangleVertices[2] - triangleVertices[0]);

		double triArea = ( e01 ^ e02).length() * 0.5; 

		MVector pv0 = triangleVertices[0] - point;
		MVector pv1 = triangleVertices[1] - point;
		MVector pv2 = triangleVertices[2] - point;

		baricentricCoords[0] = ((pv2 ^ pv1).length() * 0.5 ) / triArea;
		baricentricCoords[1] = ((pv0 ^ pv2).length() * 0.5 ) / triArea;
		baricentricCoords[2] = ((pv1 ^ pv0).length() * 0.5 ) / triArea;
	}

	MVector reflectedRay(const MVector& ligthDir, const MVector& normal )
	{
		return (ligthDir - 2 * normal * ( ligthDir * normal)).normal();
	}

	MVector halfVector(const MVector& lightDir, const MVector& viewdDir )
	{
		return (lightDir + viewdDir).normal();
	}

	MColor sumColors( const MColor& c1 , const MColor& c2 )
	{
		MColor res;
		for (int i = 0 ; i < 3; i++)
		{
			float tmp = c1[i] + c2[i];
			res[i] = (float)(tmp > 1 ? 1.0 : tmp);
		}
		return res;
	}


	void probablility()
	{
		/*
		n is number of samples so far?
		Betha - probability of error in our estimate
		sigma2 - variance of the original signal
		mu - expectation of the original signal
		E(mu n) = mu => if variance of mu n is small we are good
		var (mu n) = sigma2/n

		Condition: Prob(Var(mu n) < T) >= 1-Betha

		Stornger condition: (sigma n)2 / m(n,betha) < T
		This is because Dist[n * (signa n)2 / sigma2] is chi2


		incremental expectation
		mu n = mu (n-1) + (1/n)(last sample - mu (n-1))
		incremental variance

		Vn = (sigma n)2 * n
		Vn = V(n-1) + (sample - mu (n-1) )(sample - mu n)
		(sigma n)2 = Vn/n


		*/
	}

	double nextExpectation(double prevExpectation, int numSamples, double lastSample) {
		return prevExpectation + (1/numSamples)*(lastSample - prevExpectation);
	}

	double nextVariance(double prevVariance, double prevExpectation, double nextExpectation, int numSamples, double lastSample) {
		double prevV = prevVariance * (numSamples - 1);
		double nextV = prevV + (lastSample - prevExpectation)*(lastSample - nextExpectation);
		return nextV / numSamples;
	}

	MColor nextColorExpectation(MColor prevExpectation, int numSamples, MColor lastSample) {
		MColor nextColorE;
		for (int i = 0; i < 4; i++) {
			nextColorE[i] = nextExpectation(prevExpectation[i], numSamples, lastSample[i]);
		}
		return nextColorE;
	}

	MColor nextColorVariance(MColor prevVariance, MColor prevExpectation, MColor nextExpectation, int numSamples, MColor lastSample) {
		MColor nextColorVar;
		for (int i = 0; i < 4; i++) {
			nextColorVar[i] = nextVariance(prevVariance[i], prevExpectation[i], nextExpectation[i], numSamples, lastSample[i]); 
		}
		return nextColorVar;
	}



//scary
#pragma region 
	
	
#define X 0
#define Y 1
#define Z 2

#define CROSS(dest,v1,v2) \
	dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
	dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
	dest[2]=v1[0]*v2[1]-v1[1]*v2[0]; 

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2) \
	dest[0]=v1[0]-v2[0]; \
	dest[1]=v1[1]-v2[1]; \
	dest[2]=v1[2]-v2[2]; 

#define FINDMINMAX(x0,x1,x2,minVal,maxVal) \
	minVal = maxVal = x0;   \
	if(x1<minVal) minVal=x1;\
	if(x1>maxVal) maxVal=x1;\
	if(x2<minVal) minVal=x2;\
	if(x2>maxVal) maxVal=x2;


	/*======================== X-tests ========================*/

#define AXISTEST_X01(a, b, fa, fb, v0, v2)			   \
	p0 = a*v0[Y] - b*v0[Z];			       	   \
	p2 = a*v2[Y] - b*v2[Z];			       	   \
	if(p0<p2) {minVal=p0; maxVal=p2;} else {minVal=p2; maxVal=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(minVal>rad || maxVal<-rad) return false;


#define AXISTEST_X2(a, b, fa, fb, v0, v1)			   \
	p0 = a*v0[Y] - b*v0[Z];			           \
	p1 = a*v1[Y] - b*v1[Z];			       	   \
	if(p0<p1) {minVal=p0; maxVal=p1;} else {minVal=p1; maxVal=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(minVal>rad || maxVal<-rad) return false;



	/*======================== Y-tests ========================*/

#define AXISTEST_Y02(a, b, fa, fb, v0, v2)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p2 = -a*v2[X] + b*v2[Z];	       	       	   \
	if(p0<p2) {minVal=p0; maxVal=p2;} else {minVal=p2; maxVal=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(minVal>rad || maxVal<-rad) return false;



#define AXISTEST_Y1(a, b, fa, fb, v0, v1)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p1 = -a*v1[X] + b*v1[Z];	     	       	   \
	if(p0<p1) {minVal=p0; maxVal=p1;} else {minVal=p1; maxVal=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(minVal>rad || maxVal<-rad) return false;



	/*======================== Z-tests ========================*/



#define AXISTEST_Z12(a, b, fa, fb, v1, v2)			   \
	p1 = a*v1[X] - b*v1[Y];			           \
	p2 = a*v2[X] - b*v2[Y];			       	   \
	if(p2<p1) {minVal=p2; maxVal=p1;} else {minVal=p1; maxVal=p2;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(minVal>rad || maxVal<-rad) return false;



#define AXISTEST_Z0(a, b, fa, fb, v0, v1)			   \
	p0 = a*v0[X] - b*v0[Y];				   \
	p1 = a*v1[X] - b*v1[Y];			           \
	if(p0<p1) {minVal=p0; maxVal=p1;} else {minVal=p1; maxVal=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(minVal>rad || maxVal<-rad) return false;
	
	inline bool planeBoxOverlap(const MVector& normal,const MPoint& point, const double halfBox[3])	// -NJMP-
	{
		int q;
		double v;
		MPoint vmin, vmax;
		for(q=X;q<=Z;q++)
		{
			v=point[q];					
			if(normal[q]>0.0f)
			{
				vmin[q]=-halfBox[q] - v;	
				vmax[q]= halfBox[q] - v;	
			}
			else
			{
				vmin[q]= halfBox[q] - v;	
				vmax[q]=-halfBox[q] - v;	
			}
		}
		if((normal * vmin) > 0.0f) return false;	
		if((normal * vmax) >= 0.0f) return true;	
		return false;
	}

	bool inner_triangleBoxOverlap( const MPoint& center , const double boxhalfsize[3], const MPointArray& triangleVertices)
	{
		/*    use separating axis theorem to test overlap between triangle and box */
		/*    need to test for overlap in these directions: */
		/*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
		/*       we do not even need to test these) */
		/*    2) normal of the triangle */
		/*    3) crossproduct(edge from tri, {x,y,z}-directin) */
		/*       this gives 3x3=9 more tests */
		if(triangleVertices.length() != 3)
			return false;

		MPoint tvs[3];
		register double minVal,maxVal,p0,p1,p2,rad,fex,fey,fez;		
		MVector norm, edges[3];

		/* This is the fastest branch on Sun */
		/* move everything so that the boxcenter is in (0,0,0) */

		tvs[0] = triangleVertices[0] - center;
		tvs[1] = triangleVertices[1] - center;
		tvs[2] = triangleVertices[2] - center;

		edges[0] = tvs[1] - tvs[0];
		edges[1] = tvs[2] - tvs[1];
		edges[2] = tvs[0] - tvs[2];

		/* Bullet 3:  */
		/*  test the 9 tests first (this was faster) */

		fex = fabs(edges[0][X]);
		fey = fabs(edges[0][Y]);
		fez = fabs(edges[0][Z]);

		AXISTEST_X01(edges[0][Z], edges[0][Y], fez, fey, tvs[0], tvs[2]);
		AXISTEST_Y02(edges[0][Z], edges[0][X], fez, fex, tvs[0], tvs[2]);
		AXISTEST_Z12(edges[0][Y], edges[0][X], fey, fex, tvs[1], tvs[2]);

		fex = fabs(edges[1][X]);
		fey = fabs(edges[1][Y]);
		fez = fabs(edges[1][Z]);

		AXISTEST_X01(edges[1][Z], edges[1][Y], fez, fey, tvs[0], tvs[2]);
		AXISTEST_Y02(edges[1][Z], edges[1][X], fez, fex, tvs[0], tvs[2]);
		AXISTEST_Z0(edges[1][Y], edges[1][X], fey, fex, tvs[0], tvs[1]);

		fex = fabs(edges[2][X]);
		fey = fabs(edges[2][Y]);
		fez = fabs(edges[2][Z]);

		AXISTEST_X2(edges[2][Z], edges[2][Y], fez, fey, tvs[0] , tvs[1]);
		AXISTEST_Y1(edges[2][Z], edges[2][X], fez, fex, tvs[0], tvs[1]);
		AXISTEST_Z12(edges[2][Y], edges[2][X], fey, fex, tvs[1], tvs[2]);



		/* Bullet 1: */
		/*  first test overlap in the {x,y,z}-directions */
		/*  find minVal, maxVal of the triangle each direction, and test for overlap in */
		/*  that direction -- this is equivalent to testing a minimal AABB around */
		/*  the triangle against the AABB */

		/* test in X-direction */
		FINDMINMAX(tvs[0][X],tvs[1][X],tvs[2][X],minVal,maxVal);
		if(minVal>boxhalfsize[X] || maxVal<-boxhalfsize[X]) return false;

		/* test in Y-direction */
		FINDMINMAX(tvs[0][Y],tvs[1][Y],tvs[2][Y],minVal,maxVal);
		if(minVal>boxhalfsize[Y] || maxVal<-boxhalfsize[Y]) return false;

		/* test in Z-direction */
		FINDMINMAX(tvs[0][Z],tvs[1][Z],tvs[2][Z],minVal,maxVal);
		if(minVal>boxhalfsize[Z] || maxVal<-boxhalfsize[Z]) return false;

		/* Bullet 2: */
		/*  test if the box intersects the plane of the triangle */
		/*  compute plane equation of triangle: normal*x+d=0 */

		norm = edges[0] ^ edges[1];
		if(!planeBoxOverlap(norm, tvs[0], boxhalfsize)) return false;

		return true;   /* box and triangle overlaps */

	}


	bool triangleBoxOverlap( const MPoint& center , const double boxhalfsize[3], const MPointArray& triangleVertices)
	{
		Profiler::startTimer("SELF::triangleBoxOverlap");
		bool res = inner_triangleBoxOverlap(center, boxhalfsize, triangleVertices);
		Profiler::finishTimer("SELF::triangleBoxOverlap");
		return res;
	}

	

#pragma endregion

}

#ifndef EQNURBSNODE_H
#define EQNURBSNODE_H

#include <maya/MGlobal.h>
#include <maya/MPxNode.h>

#include <maya/MIOStream.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MArrayDataBuilder.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsCurveData.h>
#include <maya/MTime.h>
#include <maya/MPointArray.h> 
#include <maya/MDoubleArray.h>
#include <maya/MFnDoubleArrayData.h>

#include <string>
#include <math.h>

class eqNurbsNode : public MPxNode
{
public:
						eqNurbsNode();
	virtual				~eqNurbsNode(); 

	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );

	static	void*		creator();
	static	MStatus		initialize();
	void postConstructor();
 

	//static	MObject		inTime;
	static	MObject		inFFT;
	static	MObject		inFFTcompact;
	static	MObject		inFFTraw;
	static	MObject		degree;
	static	MObject		scaleX;
	static	MObject		scaleY;


	static	MObject		outCurve;
	static	MObject		outCurves;

	static	MTypeId		id;

	MObject nodeInstance;


};

#endif
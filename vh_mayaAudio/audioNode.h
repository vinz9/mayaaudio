#ifndef AUDIONODE_H
#define AUDIONODE_H

#include <maya/MPxNode.h>
#include <maya/MIOStream.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MArrayDataBuilder.h>
#include <maya/MTime.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MDoubleArray.h>


#include <string>
#include <math.h>

#include "FileRead.h"
#include "fftw3.h"

class audioNode : public MPxNode
{
public:
						audioNode();
	virtual				~audioNode(); 

	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );

	static	void*		creator();
	static	MStatus		initialize();
	//void postConstructor();
 

	static	MObject		inTime;
	static	MObject		filePath;
	static	MObject		offset;

	static	MObject		ampL;
	static	MObject		ampR;
	static	MObject		fftOutput;
	static	MObject		outFFTcompact;
	static	MObject		outFFTraw;
	static	MObject		samplesNumber;
	static	MObject		bandsNumber;
	static	MObject		ampScale;
	static	MObject     ampResponse;
	static	MObject     fftScaling;

	static	MTypeId		id;

	float lastTime;
	bool timeChanged;


	MString currentFilePath;
	//unsigned int songLength;
	long songLength;

//	bool isLoaded;

	float* fftin;
	fftwf_complex *fftout;
	fftwf_plan p;
	float* fftMag;
	MDoubleArray fftRawArray;

	double powerN;
	//int fftSamples;
	int nfftResult;
	int nBands;
	float* fftBands;
	float* fftBandsTemp;
	float curAmpScale;
	float curResponse;
	float curFftScaling;


	float halfLife;
	float sampleRate;
	float gain;

	int nBytes;
	int nChannels;
	int nSamples;
	//int chunkSize;
	int bits;

	float newAmpLeft;
	float newAmpRight;

	short* buffer;
	//int* bufferInt;

	FileRead* mySound;
	StkFrames* currentChunk;

	void peak(StkFrames* data);

	void fftw(StkFrames* buffer);
	double newton(int bands, int fftsize);

	MStatus updateFilePath(MString newFilePath);
	void updateBands(int newBandsNumber);
	void updateNSamples(int newNSamples);


};

#endif
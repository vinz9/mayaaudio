#include "audioNode.h"

MTypeId     audioNode::id( 0x80007 );

MObject     audioNode::inTime;
MObject     audioNode::filePath;
MObject     audioNode::offset;

MObject     audioNode::ampL; 
MObject     audioNode::ampR; 
MObject     audioNode::fftOutput;
MObject		audioNode::outFFTcompact;
MObject		audioNode::outFFTraw;
MObject		audioNode::samplesNumber;
MObject     audioNode::bandsNumber;
MObject     audioNode::ampScale;
MObject     audioNode::ampResponse;
MObject     audioNode::fftScaling;


//
void* audioNode::creator()
{
	return new audioNode();
}

//
MStatus audioNode::initialize()
{

	MFnUnitAttribute	uAttr;
	MFnNumericAttribute nAttr;
	MFnTypedAttribute	tAttr;
	MFnEnumAttribute	eAttr;
	MStatus				stat;

	inTime =  uAttr.create( "inTime", "t", MFnUnitAttribute::kTime, 0.0, &stat );
	CHECK_MSTATUS(stat);
	stat = uAttr.setStorable(true);


	filePath = tAttr.create( "filePath", "fp", MFnData::kString, MObject::kNullObj, &stat );
	CHECK_MSTATUS(stat);

	//offset = nAttr.create("offset", "of", MFnNumericData::kFloat, 0.0, &stat);
	offset = uAttr.create( "offset", "of", MFnUnitAttribute::kTime, 0.0, &stat );
	CHECK_MSTATUS(stat);  

	bandsNumber = nAttr.create("bandsNumber", "bn", MFnNumericData::kInt, 32, &stat);
	CHECK_MSTATUS(stat);  
	nAttr.setMin(3);
	nAttr.setMax(1024);
    
    fftScaling = eAttr.create("fftScaling", "fs", 1, &stat);
	CHECK_MSTATUS(stat);

	eAttr.addField( "No Scaling", 0 );
	eAttr.addField( "Square Root", 1 );
	eAttr.addField( "Average Log", 2 );

	samplesNumber = eAttr.create("samplesNumber", "sn", 3, &stat);
	CHECK_MSTATUS(stat);

	eAttr.addField( "512", 0 );
	eAttr.addField( "1024", 1 );
	eAttr.addField( "2048", 2 );
	eAttr.addField( "4096", 3 );
	eAttr.addField( "8192", 4 );

	ampResponse = nAttr.create("ampResponse", "ar", MFnNumericData::kFloat, 0.0, &stat);
	CHECK_MSTATUS(stat);  
	nAttr.setMin(0.0);
	nAttr.setMax(1.0);

	ampScale = nAttr.create("ampScale", "as", MFnNumericData::kFloat, 1.0, &stat);
	CHECK_MSTATUS(stat);  
	nAttr.setMin(0.01);
	nAttr.setMax(100.0);


	ampL = nAttr.create( "amplitudeLeft", "ampL", MFnNumericData::kFloat, 0.0, &stat );
	CHECK_MSTATUS(stat);
	stat = nAttr.setWritable( false );
	stat = nAttr.setStorable( false );


	ampR = nAttr.create( "amplitudeRight", "ampR", MFnNumericData::kFloat, 0.0, &stat );
	CHECK_MSTATUS(stat);
	stat = nAttr.setWritable( false );
	stat = nAttr.setStorable( false );

	fftOutput = nAttr.create( "fftOut", "fft", MFnNumericData::kFloat, 0.0, &stat );
	CHECK_MSTATUS(stat);         

	nAttr.setArray(true);
	nAttr.setUsesArrayDataBuilder( true );

	nAttr.setWritable(false);


	outFFTcompact = tAttr.create("outFFTcompact", "fftc", MFnData::kDoubleArray, MObject::kNullObj, &stat);
	CHECK_MSTATUS(stat);
	nAttr.setWritable(false);

	outFFTraw = tAttr.create("outFFTraw", "fftr", MFnData::kDoubleArray, MObject::kNullObj, &stat);
	CHECK_MSTATUS(stat);
	nAttr.setWritable(false);


	stat = addAttribute(inTime);
	stat = addAttribute(filePath);
	stat = addAttribute(bandsNumber);
	stat = addAttribute(ampScale);
	stat = addAttribute(ampResponse);
	stat = addAttribute(fftScaling);

	stat = addAttribute(ampL);
	stat = addAttribute(ampR);
	stat = addAttribute(fftOutput);
	stat = addAttribute(outFFTcompact);

	stat = addAttribute(offset);
	stat = addAttribute(outFFTraw);
	stat = addAttribute(samplesNumber);




	stat = attributeAffects(inTime, ampL);
	stat = attributeAffects(inTime, ampR);
	stat = attributeAffects(inTime, fftOutput);

	stat = attributeAffects(offset, ampL);
	stat = attributeAffects(offset, ampR);
	stat = attributeAffects(offset, fftOutput);


	stat = attributeAffects(ampScale, ampL);
	stat = attributeAffects(ampScale, ampR);
	stat = attributeAffects(ampScale, fftOutput);

	stat = attributeAffects(ampResponse, fftOutput);
	stat = attributeAffects(fftScaling, fftOutput);
	stat = attributeAffects(samplesNumber, fftOutput);

	stat = attributeAffects(samplesNumber, ampL);
	stat = attributeAffects(samplesNumber, ampR);


	stat = attributeAffects(inTime, outFFTcompact);
	stat = attributeAffects(offset, outFFTcompact);
	stat = attributeAffects(ampScale, outFFTcompact);
	stat = attributeAffects(ampResponse, outFFTcompact);
	stat = attributeAffects(fftScaling, outFFTcompact);
	stat = attributeAffects(samplesNumber, outFFTcompact);

	stat = attributeAffects(inTime, outFFTraw);
	stat = attributeAffects(offset, outFFTraw);
	//stat = attributeAffects(ampScale, outFFTraw);
	//stat = attributeAffects(ampResponse, outFFTraw);
	stat = attributeAffects(samplesNumber, outFFTraw);

	/*stat = attributeAffects(bandsNumber, ampL);
	stat = attributeAffects(bandsNumber, ampR);
	stat = attributeAffects(bandsNumber, fftOutput);
	stat = attributeAffects(bandsNumber, outFFTcompact);
	stat = attributeAffects(bandsNumber, outFFTraw);*/



	return MS::kSuccess;
} 

audioNode::audioNode() {

//	FMOD_RESULT       result;

	//std::cout << "constructor" << std::endl;

	currentFilePath = "";

	sampleRate = 44100;
	halfLife = 0.5;
	gain = 1/32768.0;
	nBands = 32;

	nBytes = 2;
	nChannels = 2;
	bits = 16;
	//nSamples = 2048;
	nSamples = 4096;

	//fftSamples = 2048;
	

	newAmpLeft = 0;
	newAmpRight = 0;

	songLength = 0;

	lastTime = -1.0;
	timeChanged = false;
	

	nfftResult = (int)(nSamples*0.5+1);

	fftRawArray.setLength(nfftResult);

	fftin = new float[nSamples];
	fftout = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nSamples);
	//p = fftwf_plan_dft_r2c_1d(nSamples, fftin, out, FFTW_ESTIMATE);
	p = fftwf_plan_dft_r2c_1d(nSamples, fftin, fftout, FFTW_MEASURE);
	fftMag = new float[nfftResult];

	fftBands = new float[nBands];
	fftBandsTemp = new float[nBands];
	powerN = log10((double)nfftResult)/nBands;

	//chunkSize = nSamples*nBytes*nChannels;
	//chunkSize = nSamples*nChannels;
	


	buffer = new short[nSamples*nChannels];

	/*result = FMOD::System_Create(&fmodSystem);
	ERRCHECK(result);

	result = fmodSystem->init(1, FMOD_INIT_NORMAL, 0);
	ERRCHECK(result);*/

}


audioNode::~audioNode() {


	//std::cout << "destructor" << std::endl;
	delete[] buffer;

	fftwf_destroy_plan(p);
	fftwf_free(fftout);
	delete[] fftin;
	delete[] fftMag;

	delete[] fftBands;
	delete[] fftBandsTemp;

	if (songLength != 0) {
		mySound->close();
		delete mySound;
		delete currentChunk;

	}


	/*result = fmodSystem->close();
	ERRCHECK(result);
	result = fmodSystem->release();
	ERRCHECK(result);*/
}


MStatus audioNode::compute (const MPlug& plug, MDataBlock& data)
{

	MStatus returnStatus;
	//FMOD_RESULT       result;

	if (plug == ampL || plug == ampR  || plug == fftOutput || plug == outFFTcompact || plug == outFFTraw) {

		MDataHandle filePathHandle = data.inputValue (filePath, &returnStatus);
		CHECK_MSTATUS( returnStatus );

		MString newFilePath = filePathHandle.asString();
		if (updateFilePath(newFilePath) != MS::kSuccess) {
			return MS::kFailure;
		}


		/*MDataHandle offsetHandle = data.inputValue (offset, &returnStatus);
		CHECK_MSTATUS( returnStatus );

		float curOffset = offsetHandle.asFloat();*/

		MDataHandle offsetHandle = data.inputValue (offset, &returnStatus);
		CHECK_MSTATUS( returnStatus );

		MTime currentOffset(offsetHandle.asTime());
		//std::cout << "time 1 " << currentOffset.value() << std::endl;
		currentOffset.setUnit(MTime::kSeconds);
		//std::cout << "time seconds " << currentOffset.value() << std::endl;

		MDataHandle nSamplesHandle = data.inputValue (samplesNumber, &returnStatus);
		CHECK_MSTATUS( returnStatus );

		int newNSamples = nSamplesHandle.asInt();
		updateNSamples(newNSamples);

		MDataHandle ampScaleHandle = data.inputValue (ampScale, &returnStatus);
		CHECK_MSTATUS( returnStatus );

		curAmpScale = ampScaleHandle.asFloat();

		MDataHandle inTimeHandle = data.inputValue (inTime, &returnStatus);
		CHECK_MSTATUS( returnStatus );

		MTime currentTime(inTimeHandle.asTime());
		currentTime.setUnit(MTime::kSeconds);

		float curTime = (float)currentTime.value();

		if (curTime != lastTime) {
			timeChanged = true;
			lastTime = curTime;
			//std::cout << "timeChanged****************" <<std::endl;
		} else {
			timeChanged = false;
		}



		//std::cout << "timeread " << currentTime << std::endl;

		//int currentSample = (int)(currentTime.value()*sampleRate-chunkSize*0.5);
		//int currentSample = (int)(currentTime.value()*sampleRate-chunkSize);
		int currentSample = (int)((currentTime.value()-currentOffset.value())*sampleRate+0.5)-(int)(nSamples*0.5);
		

		if (currentSample < 0) currentSample = 0;

		if (currentSample < (int)songLength) {


			mySound->read(*currentChunk, currentSample, 0);

			if (plug == ampL || plug == ampR) {

				newAmpLeft = 0.0;
				newAmpRight = 0.0;

				peak(currentChunk);

				newAmpLeft = newAmpLeft * curAmpScale;
				newAmpRight = newAmpRight * curAmpScale;


				MDataHandle ampLeftHandle = data.outputValue (ampL, &returnStatus);
				CHECK_MSTATUS(returnStatus);

				MDataHandle ampRightHandle = data.outputValue (ampR, &returnStatus);
				CHECK_MSTATUS(returnStatus);

				ampLeftHandle.set(newAmpLeft);
				ampRightHandle.set(newAmpRight);

			} else if (plug == fftOutput || plug == outFFTcompact || plug == outFFTraw) {
				

				MDataHandle bandsHandle = data.inputValue (bandsNumber, &returnStatus);
				CHECK_MSTATUS( returnStatus );

				int newBandsNumber = bandsHandle.asInt();
				updateBands(newBandsNumber);


				MDataHandle fftScalingHandle = data.inputValue (fftScaling, &returnStatus);
				CHECK_MSTATUS( returnStatus );

				curFftScaling = fftScalingHandle.asShort();

				MDataHandle respScaleHandle = data.inputValue (ampResponse, &returnStatus);
				CHECK_MSTATUS( returnStatus );

				curResponse = respScaleHandle.asFloat();

				if (timeChanged)
					fftw(currentChunk);


				if (plug == fftOutput) {

					MArrayDataHandle outputData = data.outputArrayValue(fftOutput, &returnStatus);

					MArrayDataBuilder outputArrayBuilder = outputData.builder(&returnStatus);
					CHECK_MSTATUS(returnStatus);

					for(int i = 0; i < nBands; i++) 
					{       
						MDataHandle outputDataElem = outputArrayBuilder.addElement(i);
						outputDataElem.set( fftBands[i] );
					}

					outputData.set(outputArrayBuilder);

				} else if (plug == outFFTcompact) {

					MDoubleArray arrayFFT;

					for(int i = 0; i < nBands; i++) 
					{       
						arrayFFT.append(fftBands[i]);
					}
					
					MFnDoubleArrayData outFFTData;
					MObject outFFT = outFFTData.create(arrayFFT, &returnStatus);

					MDataHandle outHandle = data.outputValue(outFFTcompact, &returnStatus);
					outHandle.set(outFFT);



				} else {

					//std::cout << "fftraw " << fftRawArray[0] << std::endl;
					
					MFnDoubleArrayData outFFTrawData;
					MObject outFFT = outFFTrawData.create(fftRawArray, &returnStatus);

					MDataHandle outRawHandle = data.outputValue(outFFTraw, &returnStatus);
					outRawHandle.set(outFFT);


				}

			}

		}


		data.setClean(plug);

	} else {

		return MS::kUnknownParameter;
	}

	return MS::kSuccess;

}


void audioNode::peak(StkFrames* data)
{
	float scalar = (float)pow( 0.5, 1.0/(halfLife * sampleRate));
	float inL = 0;
	float inR = 0;

	for(int i = 0; i<nSamples; i++)
	{	
		if (nChannels == 1) {

			//inL = fabs(data[i]*gain);
			inL = fabs((float)data->interpolate(i,0)*gain);
			
		if ( inL >= newAmpLeft ) {
			newAmpLeft = inL;
		} /*else {
			newAmpLeft = newAmpLeft * scalar;
			if( newAmpLeft < 0.00001 ) newAmpLeft = 0.0;
			}*/
		}

		else if (nChannels == 2) {

			//inL = fabs(data[2*i]*gain);
			inL = fabs((float)data->interpolate(i,0)*gain);
			if ( inL >= newAmpLeft ) {
				newAmpLeft = inL;
			} /*else {
			  newAmpLeft = newAmpLeft * scalar;
			  if( newAmpLeft < 0.00001 ) newAmpLeft = 0.0;
			  }*/


			//inR = fabs(data[2*i+1]*gain);
			inL = fabs((float)data->interpolate(i,1)*gain);
			if ( inR >= newAmpRight )
			{
				newAmpRight = inR;
			} /*else {
			  newAmpRight = newAmpRight * scalar;
			  if( newAmpRight < 0.00001 ) newAmpRight = 0.0;
			  }*/
		} else {
			std::cout << "nChannels error" <<std::endl;
		}
	}
}



/*void audioNode::ERRCHECK(FMOD_RESULT result)
{
	if (result != FMOD_OK)
	{
		std::cout << "FMOD error!" << FMOD_ErrorString(result) << std::endl;

	}
}*/

MStatus audioNode::updateFilePath(MString newFilePath) {

	//FMOD_RESULT result;
	
	if (newFilePath != currentFilePath) {

		//std::cout << "newfilepath" << newFilePath << std::endl;

		if (songLength != 0) {

			mySound->close();
			delete mySound;
			delete currentChunk;
		}

		//result = fmodSystem->createStream(newFilePath.asChar(), FMOD_OPENONLY | FMOD_ACCURATETIME, 0, &sound);
		mySound = new FileRead(newFilePath.asChar());
		songLength = mySound->fileSize();
		
		nChannels = mySound->channels();
		currentChunk = new StkFrames(nSamples,nChannels);
	

	} else if (currentFilePath == "") {
			std::cout << "No File" << std::endl;
			return MS::kFailure;
	}

	return MS::kSuccess;
}

void audioNode::updateBands(int newBandsNumber){

	if (newBandsNumber != nBands) {

		nBands = newBandsNumber;
		powerN = log10((double)nfftResult)/nBands;

		delete[] fftBands;
		delete[] fftBandsTemp;
		fftBands = new float[nBands];
		fftBandsTemp = new float[nBands];
	}
}

void audioNode::updateNSamples(int newNSamples){

	int tempSamples = (int)(pow((double)2, 9+newNSamples));
	if (tempSamples != nSamples) {
		nSamples = tempSamples;

		delete[] buffer;

		fftwf_destroy_plan(p);
		fftwf_free(fftout);
		delete[] fftin;
		delete[] fftMag;

		delete[] fftBands;
		delete[] fftBandsTemp;

		nfftResult = (int)(nSamples*0.5+1);
		fftRawArray.setLength(nfftResult);

		fftin = new float[nSamples];
		fftout = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * nSamples);
		//p = fftwf_plan_dft_r2c_1d(nSamples, fftin, out, FFTW_ESTIMATE);
		p = fftwf_plan_dft_r2c_1d(nSamples, fftin, fftout, FFTW_MEASURE);
		fftMag = new float[nfftResult];

		fftBands = new float[nBands];
		fftBandsTemp = new float[nBands];
		powerN = log10((double)nfftResult)/nBands;

//		chunkSize = nSamples*nBytes*nChannels;

		buffer = new short[nSamples*nChannels];


	}


}

void audioNode::fftw(StkFrames* buffer)
{

	//std::cout << "fftwread" << std::endl;

	double window, imag, real;

	for (int k = 0; k < nSamples; k++) {
		window = -.5f*cos(2.*M_PI*k/nSamples)+.5f;
		if (nChannels == 2) {
		fftin[k] = (float)(((float)buffer->interpolate(k,0)+(float)buffer->interpolate(k,1))*0.5*gain * window);
		//fftin[k] = (float)(buffer[2*k]*0.5*gain * window);

		} else if (nChannels == 1) {
		fftin[k] = (float)((float)buffer->interpolate(k,0)*gain * window);
		}
	}

	fftwf_execute(p);

	for (int k = 0; k < nfftResult; k++) {

		real = fftout[k][0];
		imag = fftout[k][1];

		fftMag[k] = (float)sqrt(real*real + imag*imag);
		//fftRawArray.append(fftMag[k]*ampScale);
		fftRawArray.set(fftMag[k]*curAmpScale, k);
	}


	int freq = 0;
	int freqlast = 0;
	float result = 0;

	for(int i=1; i<=nBands; i++)
	{
		float sum = 0.0f;
		float sumAvg = 0.0f;
		
		freqlast = freq;
		freq = (int)(pow((double)10, i*powerN)+0.5);
		if (freq >= nfftResult) freq = nfftResult - 1;

		if (freqlast == freq) freqlast = freq-1;
		int bsize = freq-freqlast;

		for(;freqlast < freq; freqlast++){
			sum += fftMag[freqlast];
		}


		if (curFftScaling == 0) {
			sum = (float)(sum/500.0);
			result = sum;
		}

		else if (curFftScaling == 1) {
			
			//sum = (float)(sum/500.0);
			sum = (float)(sum/500.0);
			//std::cout << "sum " << sum << std::endl;
			//sum = sum/32.0*nBands;
			//std::cout << "sumnorm " << sum << std::endl;
			result = (float)sqrt((double)sum);
		}

		else if (curFftScaling == 2) {

			sumAvg = (float)sum/bsize;

			if (sumAvg > 0.01)
				sumAvg = (float)(log((double)sumAvg) + 4.6);
			else
				sumAvg = 0;

			result = (float)(sumAvg/10.0);
		}


		//result = result*curAmpScale;

		if (curResponse == 0)
			fftBandsTemp[i-1] = result;
		else {
			if ( result > fftBandsTemp[i-1] || (fftBandsTemp[i-1]  - result < 0.001) )
				fftBandsTemp[i-1] = result;
			else
			{
				fftBandsTemp[i-1] = fftBandsTemp[i-1] * curResponse;
				if( fftBandsTemp[i-1] < 0.00001 ) fftBandsTemp[i-1] = 0.0;
			}
		}

		fftBands[i-1] = fftBandsTemp[i-1]*curAmpScale;
	}


}

double audioNode::newton(int nBands, int fftsize)
{
	//cout << "start newton" << endl;
	double x = 0.0; //power?
	int b = nBands; //number of output nBands
	int fs = fftsize; //number of fft values
	//cout << fftsize << " fft-values into " << nBands << " nBands" << endl;
	double ic= 0.0; //substitute for int i
	double e; //error -> runs to 0
	double fx; //f(x)
	double fx1; //f'(x)
	double p, p2;
	double l;

	fx = -fs; //f(x)
	fx1 = 0.0; //f'(x)

	//values to get started:
	x = 1.0;
	e = 1.0;
	p = 1.0;
	l = 1.0;
	p2 = pow(2.0,(-15.0));

	int count = 0;

	while(abs(e)> p2)                       //Do Until Abs(dh) < 2 ^ (-15)
	{
		count++;
		int band = 0;
		for(int i = 1; i<= b; i++)
		{
			p = pow((double)i,x);
			l = log((double)i);
			fx = fx + p;
			fx1 = fx1 + p*l;
		}

		e = -fx/fx1;
		x += e;
		fx = -fs;
		fx1 = 0;
	}

	return x;
}
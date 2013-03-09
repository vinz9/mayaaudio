#include "eqNurbsNode.h"

MTypeId     eqNurbsNode::id( 0x80011 );

//MObject     eqNurbsNode::inTime;
MObject     eqNurbsNode::inFFT;
MObject     eqNurbsNode::inFFTcompact;
MObject     eqNurbsNode::inFFTraw;
MObject     eqNurbsNode::degree;
MObject     eqNurbsNode::scaleX;
MObject     eqNurbsNode::scaleY;

MObject     eqNurbsNode::outCurve; 
MObject     eqNurbsNode::outCurves; 

//
void* eqNurbsNode::creator()
{
	return new eqNurbsNode();
}

//
MStatus eqNurbsNode::initialize()
{

	MFnUnitAttribute	uAttr;
	MFnNumericAttribute nAttr;
	MFnTypedAttribute	tAttr;
	MFnEnumAttribute	eAttr;
	MStatus				stat;

	/*inTime =  uAttr.create( "inTime", "t", MFnUnitAttribute::kTime, 0.0, &stat );
	CHECK_MSTATUS(stat);
	stat = uAttr.setStorable(true);*/

	inFFT = nAttr.create( "inFFT", "ifft", MFnNumericData::kFloat, 0.0, &stat );
	CHECK_MSTATUS(stat);         

	nAttr.setArray(true);
	//stat = nAttr.setKeyable(true);
	CHECK_MSTATUS(stat);

	inFFTcompact = tAttr.create( "inFFTcompact", "ifftc", MFnData::kDoubleArray, MObject::kNullObj, &stat );
	CHECK_MSTATUS(stat); 

	inFFTraw = tAttr.create( "inFFTraw", "ifftr", MFnData::kDoubleArray, MObject::kNullObj, &stat );
	CHECK_MSTATUS(stat); 

	//stat = tAttr.setKeyable(true);
	//CHECK_MSTATUS(stat);

	//nAttr.setUsesArrayDataBuilder( true );
	//nAttr.setWritable(false);

	degree = nAttr.create( "degree", "deg", MFnNumericData::kInt, 2, &stat );
	CHECK_MSTATUS(stat); 
	nAttr.setMin(1);
	nAttr.setMax(3);

	scaleX = nAttr.create( "scaleX", "sx", MFnNumericData::kFloat, 1.0, &stat );
	CHECK_MSTATUS(stat); 
	nAttr.setMin(0.01);
	nAttr.setMax(20.0);

	scaleY = nAttr.create( "scaleY", "sy", MFnNumericData::kFloat, 0.50, &stat );
	CHECK_MSTATUS(stat); 
	nAttr.setMin(0.01);
	nAttr.setMax(20.0);

	outCurve = tAttr.create("outCurve", "oc", MFnNurbsCurveData::kNurbsCurve, &stat );
	CHECK_MSTATUS(stat);

	tAttr.setWritable(false);

		outCurves = tAttr.create("outCurves", "ocs", MFnNurbsCurveData::kNurbsCurve, &stat );
	CHECK_MSTATUS(stat);

	tAttr.setArray(true);
	tAttr.setUsesArrayDataBuilder(true);

	//stat = addAttribute(inTime);
	stat = addAttribute(inFFT);
	stat = addAttribute(inFFTcompact);
	stat = addAttribute(inFFTraw);
	stat = addAttribute(degree);
	stat = addAttribute(scaleX);
	stat = addAttribute(scaleY);
	stat = addAttribute(outCurve);	

	stat = attributeAffects(inFFT, outCurve);
	stat = attributeAffects(inFFTcompact, outCurve);
	stat = attributeAffects(inFFTraw, outCurve);
	stat = attributeAffects(degree, outCurve);
	stat = attributeAffects(scaleX, outCurve);
	stat = attributeAffects(scaleY, outCurve);

	stat = attributeAffects(inFFTcompact, outCurves);
	stat = attributeAffects(degree, outCurves);
	stat = attributeAffects(scaleX, outCurves);
	stat = attributeAffects(scaleY, outCurves);

	return MS::kSuccess;
} 



eqNurbsNode::eqNurbsNode() {

}

void eqNurbsNode::postConstructor() {

	nodeInstance = thisMObject();

}

eqNurbsNode::~eqNurbsNode() {


}


MStatus eqNurbsNode::compute (const MPlug& plug, MDataBlock& data)
{

	MStatus returnStatus;

	if (plug == outCurve || plug == outCurves) {

		MPlug plugFFTcompact(nodeInstance, inFFTcompact);
		MPlug plugFFTraw(nodeInstance, inFFTraw);

		MPointArray cvArray;
		MDoubleArray knotArray;

		MDataHandle degHandle = data.inputValue(degree, &returnStatus);
		CHECK_MSTATUS(returnStatus);

		int curDeg = degHandle.asInt();

		MDataHandle scaleXHandle = data.inputValue(scaleX, &returnStatus);
		CHECK_MSTATUS(returnStatus);

		float curScaleX = scaleXHandle.asFloat();

		MDataHandle scaleYHandle = data.inputValue(scaleY, &returnStatus);
		CHECK_MSTATUS(returnStatus);

		float curScaleY = scaleYHandle.asFloat();



		if (plug == outCurve) {


			if (plugFFTcompact.isConnected()) {

				//std::cout << "fftcompact" << std::endl;

				MDataHandle inFFTcompactHandle = data.inputValue(inFFTcompact, &returnStatus);
				CHECK_MSTATUS(returnStatus);

				MFnDoubleArrayData fftData(inFFTcompactHandle.data(), &returnStatus);
				CHECK_MSTATUS(returnStatus);

				int nBands = fftData.length();

				for(int i = 0; i < nBands; i++) 
				{       
					cvArray.append(fftData[i]*curScaleX,(float)i/nBands*32.0*curScaleY,0);
					//std::cout << "fftdata" << fftData[i] << std::endl;

					if (i < nBands - (curDeg - 1)) {
						knotArray.append(i);

						if (i == 0 || i == (nBands-curDeg)) {
							for (int j = 0; j<curDeg-1; j++) {
								knotArray.append(i);
							}
						}
					}
				}

				//std::cout << "knotarray" << knotArray.length() << std::endl;

			} else 	if (plugFFTraw.isConnected()) {

				//std::cout << "fftraw" << std::endl;

				MDataHandle inFFTrawHandle = data.inputValue(inFFTraw, &returnStatus);
				CHECK_MSTATUS(returnStatus);

				MFnDoubleArrayData fftRawData(inFFTrawHandle.data(), &returnStatus);
				CHECK_MSTATUS(returnStatus);

				int fftSize = fftRawData.length();

				//std::cout << "fftsize " << fftSize << std::endl;

				for(int i = 0; i < fftSize; i++) 
				{       

					double x = sqrt((double)fftRawData[i])*pow((double)i,0.25)*curScaleX*0.5;

					//double t1 = log10((double)(i+1));
					//double t2 = log10((double)(fftSize));

					double y = ( log10((double)(1+i)) / log10((double)(fftSize)))*curScaleY*10.0;

					//double y = t1/t2 * curScaleY;

					//std::cout << "i :" << i << " x :" << x << " y :" << y << std::endl;

					cvArray.append(x,y,0);


					if (i < fftSize - (curDeg - 1)) {
						knotArray.append(i);

						if (i == 0 || i == (fftSize-curDeg)) {
							for (int j = 0; j<curDeg-1; j++) {
								knotArray.append(i);
							}
						}
					}
				}

			} else {

				//std::cout << "fftarray" << std::endl;

				MArrayDataHandle inFFTArray = data.inputArrayValue(inFFT, &returnStatus);
				CHECK_MSTATUS(returnStatus);

				//MArrayDataBuilder inFFTArrayBuilder = inFFTArray.builder(&returnStatus);
				//CHECK_MSTATUS(returnStatus);

				//int nBands = inFFTArrayBuilder.elementCount(&returnStatus);

				int nBands = inFFTArray.elementCount(&returnStatus);

				inFFTArray.jumpToElement(0);

				for(int i = 0; i < nBands; i++) 
				{       
					MDataHandle fftBandHandle = inFFTArray.inputValue(&returnStatus);
					cvArray.append(fftBandHandle.asFloat()*curScaleX,(float)i/nBands*32.0*curScaleY,0);

					if (i < nBands - (curDeg - 1)) {
						knotArray.append(i);

						for (int j = 0; j<curDeg; j++) {
							knotArray.append(i);
						}
					}
					inFFTArray.next();
				}

			}


			MFnNurbsCurveData dataCreator;
			MObject outCurveData = dataCreator.create();

			MFnNurbsCurve curveFn;

			MObject curve = curveFn.create(cvArray, knotArray, curDeg, 
				MFnNurbsCurve::kOpen, false, false, outCurveData, &returnStatus); 
			CHECK_MSTATUS(returnStatus);


			MDataHandle outCurveHandle = data.outputValue(outCurve, &returnStatus);
			CHECK_MSTATUS(returnStatus);

			outCurveHandle.set(outCurveData);



			data.setClean(plug);
		} else {

			std::cout << "loft " << std::endl;

			MDataHandle inFFTcompactHandle = data.inputValue(inFFTcompact, &returnStatus);
			CHECK_MSTATUS(returnStatus);

			MFnDoubleArrayData fftData(inFFTcompactHandle.data(), &returnStatus);
			CHECK_MSTATUS(returnStatus);

			int nBands = fftData.length();

			for(int i = 0; i < nBands; i++) 
			{       
				cvArray.append(fftData[i]*curScaleX,(float)i/nBands*32.0*curScaleY,0);
				//std::cout << "fftdata" << fftData[i] << std::endl;

				if (i < nBands - (curDeg - 1)) {
					knotArray.append(i);

					if (i == 0 || i == (nBands-curDeg)) {
						for (int j = 0; j<curDeg-1; j++) {
							knotArray.append(i);
						}
					}
				}
			}

			MFnNurbsCurveData dataCreator;
			MObject outCurveData = dataCreator.create();

			MFnNurbsCurve curveFn;

			MObject curve = curveFn.create(cvArray, knotArray, curDeg, 
				MFnNurbsCurve::kOpen, false, false, outCurveData, &returnStatus); 
			CHECK_MSTATUS(returnStatus);


			MDataHandle outCurveHandle = data.outputValue(outCurve, &returnStatus);
			CHECK_MSTATUS(returnStatus);

			outCurveHandle.set(outCurveData);



			data.setClean(plug);



			}

	} else {
		return MS::kUnknownParameter;
	}

	return MS::kSuccess;
}
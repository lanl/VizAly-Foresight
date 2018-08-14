/**
 *  @file sz_float.c
 *  @author Sheng Di and Dingwen Tao
 *  @date Aug, 2016
 *  @brief SZ_Init, Compression and Decompression functions
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "sz.h"
#include "CompressElement.h"
#include "DynamicByteArray.h"
#include "DynamicIntArray.h"
#include "TightDataPointStorageF.h"
#include "zlib.h"
#include "rw.h"

unsigned int optimize_intervals_float_1D(float *oriData, int dataLength, double realPrecision)
{	
	int i = 0;
	unsigned long radiusIndex;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	int totalSampleSize = dataLength/sampleDistance;
	for(i=2;i<dataLength;i++)
	{
		if(i%sampleDistance==0)
		{
			//pred_value = 2*oriData[i-1] - oriData[i-2];
			pred_value = oriData[i-1];
			pred_err = fabs(pred_value - oriData[i]);
			radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
			if(radiusIndex>=maxRangeRadius)
				radiusIndex = maxRangeRadius - 1;			
			intervals[radiusIndex]++;
		}
	}
	//compute the appropriate number
	int targetCount = (int)(totalSampleSize*predThreshold);
	int sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
		
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);
	
	if(powerOf2<32)
		powerOf2 = 32;
	
	free(intervals);
	//printf("accIntervals=%d, powerOf2=%d\n", accIntervals, powerOf2);
	return powerOf2;
}

unsigned int optimize_intervals_float_2D(float *oriData, int r1, int r2, double realPrecision)
{	
	int i,j, index;
	unsigned long radiusIndex;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	int totalSampleSize = r1*r2/sampleDistance;
	for(i=1;i<r1;i++)
	{
		for(j=1;j<r2;j++)
		{
			if((i+j)%sampleDistance==0)
			{
				index = i*r2+j;
				pred_value = oriData[index-1] + oriData[index-r2] - oriData[index-r2-1];
				pred_err = fabs(pred_value - oriData[index]);
				radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
				if(radiusIndex>=maxRangeRadius)
					radiusIndex = maxRangeRadius - 1;
				intervals[radiusIndex]++;
			}			
		}
	}
	//compute the appropriate number
	int targetCount = (int)(totalSampleSize*predThreshold);
	int sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	//printf("maxRangeRadius = %d, accIntervals=%d, powerOf2=%d\n", maxRangeRadius, accIntervals, powerOf2);
	return powerOf2;
}

unsigned int optimize_intervals_float_3D(float *oriData, int r1, int r2, int r3, double realPrecision)
{	
	int i,j,k, index;
	unsigned long radiusIndex;
	int r23=r2*r3;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	int totalSampleSize = r1*r2*r3/sampleDistance;
	for(i=1;i<r1;i++)
	{
		for(j=1;j<r2;j++)
		{
			for(k=1;k<r3;k++)
			{			
				if((i+j+k)%sampleDistance==0)
				{
					index = i*r23+j*r3+k;
					pred_value = oriData[index-1] + oriData[index-r3] + oriData[index-r23] 
					- oriData[index-1-r23] - oriData[index-r3-1] - oriData[index-r3-r23] + oriData[index-r3-r23-1];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
					if(radiusIndex>=maxRangeRadius)
						radiusIndex = maxRangeRadius - 1;
					intervals[radiusIndex]++;
				}
			}
		}
	}
	//compute the appropriate number
	int targetCount = (int)(totalSampleSize*predThreshold);
	int sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;
	
	free(intervals);
	//printf("accIntervals=%d, powerOf2=%d\n", accIntervals, powerOf2);
	return powerOf2;
}


unsigned int optimize_intervals_float_4D(float *oriData, int r1, int r2, int r3, int r4, double realPrecision)
{
	int i,j,k,l, index;
	unsigned long radiusIndex;
	int r234=r2*r3*r4;
	int r34=r3*r4;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	int totalSampleSize = r1*r2*r3*r4/sampleDistance;
	for(i=1;i<r1;i++)
	{
		for(j=1;j<r2;j++)
		{
			for(k=1;k<r3;k++)
			{
				for (l=1;l<r4;l++)
				{
					if((i+j+k+l)%sampleDistance==0)
					{
						index = i*r234+j*r34+k*r4+l;
						pred_value = oriData[index-1] + oriData[index-r3] + oriData[index-r34]
								- oriData[index-1-r34] - oriData[index-r4-1] - oriData[index-r4-r34] + oriData[index-r4-r34-1];
						pred_err = fabs(pred_value - oriData[index]);
						radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
						if(radiusIndex>=maxRangeRadius)
							radiusIndex = maxRangeRadius - 1;
						intervals[radiusIndex]++;
					}
				}
			}
		}
	}
	//compute the appropriate number
	int targetCount = (int)(totalSampleSize*predThreshold);
	int sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	return powerOf2;
}

TightDataPointStorageF* SZ_compress_float_1D_MDQ(float *oriData, 
int dataLength, double realPrecision, float valueRangeSize, float medianValue_f)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
		quantization_intervals = optimize_intervals_float_1D(oriData, dataLength, realPrecision);
	else
		quantization_intervals = intvCapacity;
	updateQuantizationInfo(quantization_intervals);	
	//clearHuffmanMem();
	int i, reqLength;
	float medianValue = medianValue_f;
	short reqExpo = getPrecisionReqLength_float((float)realPrecision);
	short radExpo = getExponent_float(valueRangeSize/2);
	
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);	

	int* type = (int*) malloc(dataLength*sizeof(int));
		
	float* spaceFillingValue = oriData; //
	
	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);
	
	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);
	
	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);
	
	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);
	
	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);
	
	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;
	float last3CmprsData[3] = {0};

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
				
	//add the first data	
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_float(last3CmprsData, vce->data);
	//printf("%.30G\n",last3CmprsData[0]);	
		
	//add the second data
	type[1] = 0;
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);			
	compressSingleFloatValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_float(last3CmprsData, vce->data);
	//printf("%.30G\n",last3CmprsData[0]);	
	
	int state;
	float lcf, qcf;	
	double checkRadius;
	float curData;
	float pred;
	float predAbsErr;
	float min_pred, minErr, minIndex;
	int a = 0;		
	checkRadius = (intvCapacity-1)*realPrecision;
	double interval = 2*realPrecision;
	
	for(i=2;i<dataLength;i++)
	{
//		if(i==6)
//			printf("i=%d\n", i);
		curData = spaceFillingValue[i];
		//pred = 2*last3CmprsData[0] - last3CmprsData[1];
		pred = last3CmprsData[0];
		predAbsErr = fabs(curData - pred);	
		if(predAbsErr<=checkRadius)
		{
			state = (predAbsErr/realPrecision+1)/2;
			if(curData>=pred)
			{
				type[i] = intvRadius+state;
				pred = pred + state*interval;
			}
			else //curData<pred
			{
				type[i] = intvRadius-state;
				pred = pred - state*interval;
			}
/*			if(type[i]==0)
				printf("err:type[%d]=0\n", i);*/
			listAdd_float(last3CmprsData, pred);					
			continue;
		}
		
		//unpredictable data processing		
		type[i] = 0;
		addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
		
		compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);

		listAdd_float(last3CmprsData, vce->data);	
	}//end of for
		
//	char* expSegmentsInBytes;
//	int expSegmentsInBytes_size = convertESCToBytes(esc, &expSegmentsInBytes);
	int exactDataNum = exactLeadNumArray->size;
	
	TightDataPointStorageF* tdps;
			
	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum, 
			type, exactMidByteArray->array, exactMidByteArray->size,  
			exactLeadNumArray->array,  
			resiBitArray->array, resiBitArray->size, 
			resiBitLengthArray->array, resiBitLengthArray->size, 
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

//sdi:Debug
/*	int sum =0;
	for(i=0;i<dataLength;i++)
		if(type[i]==0) sum++;
	printf("opt_quantizations=%d, exactDataNum=%d, sum=%d\n",quantization_intervals, exactDataNum, sum);*/

//	writeShortData(type, dataLength, "compressStateBytes.sb");
//	unsigned short type_[dataLength];
//	SZ_Reset();
//	decode_withTree(tdps->typeArray, tdps->typeArray_size, type_);	
//	printf("tdps->typeArray_size=%d\n", tdps->typeArray_size);
	
//	printf("exactDataNum=%d, expSegmentsInBytes_size=%d, exactMidByteArray->size=%d,resiBitLengthArray->size=%d\n", 
//			exactDataNum, expSegmentsInBytes_size, exactMidByteArray->size, resiBitLengthArray->size);
	
//	for(i = 3800;i<3844;i++)
//		printf("exactLeadNumArray->array[%d]=%d\n",i,exactLeadNumArray->array[i]);
	
	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);	
	free(vce);
	free(lce);	
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);
	
	return tdps;
}

void SZ_compress_args_float_StoreOriData(float* oriData, int dataLength, TightDataPointStorageF* tdps, 
unsigned char** newByteData, int *outSize)
{
	int floatSize=sizeof(float);	
	int k = 0, i;
	tdps->isLossless = 1;
	int totalByteLength = 3 + 4 + 1 + floatSize*dataLength;
	*newByteData = (unsigned char*)malloc(totalByteLength);
	
	unsigned char dsLengthBytes[4];
	intToBytes_bigEndian(dsLengthBytes, dataLength);//4
	for (i = 0; i < 3; i++)//3
		(*newByteData)[k++] = versionNumber[i];
	for (i = 0; i < 4; i++)//4
		(*newByteData)[k++] = dsLengthBytes[i];
	(*newByteData)[k++] = 16;	//=00010000	
	
	if(sysEndianType==BIG_ENDIAN_SYSTEM)
		memcpy((*newByteData)+8, oriData, dataLength*floatSize);
	else
	{
		unsigned char* p = (*newByteData)+8;
		for(i=0;i<dataLength;i++,p+=floatSize)
			floatToBytes(p, oriData[i]);
	}
	*outSize = totalByteLength;
}

void SZ_compress_args_float_NoCkRngeNoGzip_1D(unsigned char** newByteData, float *oriData, 
int dataLength, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f)
{
	SZ_Reset();	
	TightDataPointStorageF* tdps = SZ_compress_float_1D_MDQ(oriData, dataLength, realPrecision, valueRangeSize, medianValue_f);
			
	//TODO: return bytes....
	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);
	
	if(*outSize>dataLength*sizeof(float))
		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);
	
	free_TightDataPointStorageF(tdps);
}

TightDataPointStorageF* SZ_compress_float_2D_MDQ(float *oriData, int r1, int r2, double realPrecision, float valueRangeSize, float medianValue_f)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_2D(oriData, r1, r2, realPrecision);
		updateQuantizationInfo(quantization_intervals);
	}	
	else
		quantization_intervals = intvCapacity;
	
	int i,j, reqLength;
	float pred1D, pred2D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;
		
	int dataLength = r1*r2;	
	
	P0 = (float*)malloc(r2*sizeof(float));
	memset(P0, 0, r2*sizeof(float));
	P1 = (float*)malloc(r2*sizeof(float));
	memset(P1, 0, r2*sizeof(float));
		
	float medianValue = medianValue_f;
	short radExpo = getExponent_float(valueRangeSize/2);
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);	

	int* type = (int*) malloc(dataLength*sizeof(int));
	//type[dataLength]=0;
		
	float* spaceFillingValue = oriData; //
	
	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);
	
	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);
	
	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);
	
	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);
	
	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);
	
	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
			
	/* Process Row-0 data 0*/
	type[0] = 0;
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[0] = vce->data;

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	diff = spaceFillingValue[1] - pred1D;

	itvNum =  fabs(diff)/realPrecision + 1;

	if (itvNum < intvCapacity)
	{
		if (diff < 0) itvNum = -itvNum;
		type[1] = (int) (itvNum/2) + intvRadius;
		P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;
	}
	else
	{
		type[1] = 0;

		addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
		compressSingleFloatValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[1] = vce->data;
	}

    /* Process Row-0 data 2 --> data r2-1 */
	for (j = 2; j < r2; j++)
	{
		pred1D = 2*P1[j-1] - P1[j-2];
		diff = spaceFillingValue[j] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[j] = (int) (itvNum/2) + intvRadius;
			P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
		}
		else
		{
			type[j] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[j], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[j] = vce->data;
		}
	}

	/* Process Row-1 --> Row-r1-1 */
	int index;
	for (i = 1; i < r1; i++)
	{	
		/* Process row-i data 0 */
		index = i*r2;
		pred1D = P1[0];
		diff = spaceFillingValue[index] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P0[0] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
		}
		else
		{
			type[index] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[0] = vce->data;
		}
									
		/* Process row-i data 1 --> r2-1*/
		for (j = 1; j < r2; j++)
		{
			index = i*r2+j;		

			pred2D = P0[j-1] + P1[j] - P1[j-1];

			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				type[index] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[j] = vce->data;
			}
		}

		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}
	
	if(r2!=1)
		free(P0);
	free(P1);			
	int exactDataNum = exactLeadNumArray->size;
	
	TightDataPointStorageF* tdps;
			
	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum, 
			type, exactMidByteArray->array, exactMidByteArray->size,  
			exactLeadNumArray->array,  
			resiBitArray->array, resiBitArray->size, 
			resiBitLengthArray->array, resiBitLengthArray->size, 
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

//	printf("exactDataNum=%d, expSegmentsInBytes_size=%d, exactMidByteArray->size=%d,resiBitLengthArray->size=%d\n", 
//			exactDataNum, expSegmentsInBytes_size, exactMidByteArray->size, resiBitLengthArray->size);
	
//	for(i = 3800;i<3844;i++)
//		printf("exactLeadNumArray->array[%d]=%d\n",i,exactLeadNumArray->array[i]);
	
	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);
	
	return tdps;	
}

/**
 * 
 * Note: @r1 is high dimension
 * 		 @r2 is low dimension 
 * */
void SZ_compress_args_float_NoCkRngeNoGzip_2D(unsigned char** newByteData, float *oriData, int r1, int r2, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f)
{
	SZ_Reset();	
		
	TightDataPointStorageF* tdps = SZ_compress_float_2D_MDQ(oriData, r1, r2, realPrecision, valueRangeSize, medianValue_f);

	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);

	int dataLength = r1*r2;
	if(*outSize>dataLength*sizeof(float))
		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);
	
	free_TightDataPointStorageF(tdps);	
}

TightDataPointStorageF* SZ_compress_float_3D_MDQ(float *oriData, int r1, int r2, int r3, double realPrecision, float valueRangeSize, float medianValue_f)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_3D(oriData, r1, r2, r3, realPrecision);
		updateQuantizationInfo(quantization_intervals);
	}	
	else
		quantization_intervals = intvCapacity;
	//clearHuffmanMem();
	int i,j,k, reqLength;
	float pred1D, pred2D, pred3D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;

	int dataLength = r1*r2*r3;
	int r23 = r2*r3;
	P0 = (float*)malloc(r23*sizeof(float));
	P1 = (float*)malloc(r23*sizeof(float));

	float medianValue = medianValue_f;
	short radExpo = getExponent_float(valueRangeSize/2);
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);	

	int* type = (int*) malloc(dataLength*sizeof(int));

	float* spaceFillingValue = oriData; //
	
	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);
	
	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));


	///////////////////////////	Process layer-0 ///////////////////////////
	/* Process Row-0 data 0*/
	type[0] = 0;
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[0] = vce->data;

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	diff = spaceFillingValue[1] - pred1D;

	itvNum = fabs(diff)/realPrecision + 1;

	if (itvNum < intvCapacity)
	{
		if (diff < 0) itvNum = -itvNum;
		type[1] = (int) (itvNum/2) + intvRadius;
		P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;
	}
	else
	{
		type[1] = 0;

		addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
		compressSingleFloatValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[1] = vce->data;
	}

    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++)
	{
		pred1D = 2*P1[j-1] - P1[j-2];
		diff = spaceFillingValue[j] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[j] = (int) (itvNum/2) + intvRadius;
			P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
		}
		else
		{
			type[j] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[j], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[j] = vce->data;
		}
	}

	/* Process Row-1 --> Row-r2-1 */
	int index;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	
		pred1D = P1[index-r3];
		diff = spaceFillingValue[index] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P1[index] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
		}
		else
		{
			type[index] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index] = vce->data;
		}

		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++)
		{
			index = i*r3+j;			
			pred2D = P1[index-1] + P1[index-r3] - P1[index-r3-1];

			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				type[index] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index] = vce->data;
			}
		}
	}


	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		/* Process Row-0 data 0*/
		index = k*r23;		
		pred1D = P1[0];
		diff = spaceFillingValue[index] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P0[0] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
		}
		else
		{
			type[index] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[0] = vce->data;
		}

	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			//index = k*r2*r3+j;		
			index ++;
			pred2D = P0[j-1] + P1[j] - P1[j-1];
			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
/*				if(type[index]==0)
					printf("err:type[%d]=0, index4\n", index);					*/
			}
			else
			{
				type[index] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[j] = vce->data;
			}
		}

	    /* Process Row-1 --> Row-r2-1 */
		int index2D;
		for (i = 1; i < r2; i++)
		{
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			index2D = i*r3;		
			pred2D = P0[index2D-r3] + P1[index2D] - P1[index2D-r3];
			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				type[index] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[index2D] = vce->data;
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
//				if(k==63&&i==43&&j==27)
//					printf("i=%d\n", i);
				//index = k*r2*r3 + i*r3 + j;			
				index ++;
				index2D = i*r3 + j;
				pred3D = P0[index2D-1] + P0[index2D-r3]+ P1[index2D] - P0[index2D-r3-1] - P1[index2D-r3] - P1[index2D-1] + P1[index2D-r3-1];
				diff = spaceFillingValue[index] - pred3D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
				}
				else
				{
					type[index] = 0;

					addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
					compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}
			}
		}

		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}
	if(r23!=1)
		free(P0);
	free(P1);
	int exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitLengthArray->array, resiBitLengthArray->size, 
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

//sdi:Debug
/*	int sum =0;
	for(i=0;i<dataLength;i++)
		if(type[i]==0) sum++;
	printf("opt_quantizations=%d, exactDataNum=%d, sum=%d\n",quantization_intervals, exactDataNum, sum);
*/

//	printf("exactDataNum=%d, expSegmentsInBytes_size=%d, exactMidByteArray->size=%d,resiBitLengthArray->size=%d\n",
//			exactDataNum, expSegmentsInBytes_size, exactMidByteArray->size, resiBitLengthArray->size);

	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);	
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);
	
	return tdps;	
}


void SZ_compress_args_float_NoCkRngeNoGzip_3D(unsigned char** newByteData, float *oriData, int r1, int r2, int r3, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f)
{
	SZ_Reset();	

	TightDataPointStorageF* tdps = SZ_compress_float_3D_MDQ(oriData, r1, r2, r3, realPrecision, valueRangeSize, medianValue_f);

	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);

	int dataLength = r1*r2*r3;
	if(*outSize>dataLength*sizeof(float))
		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);

	free_TightDataPointStorageF(tdps);
}


TightDataPointStorageF* SZ_compress_float_4D_MDQ(float *oriData, int r1, int r2, int r3, int r4, double realPrecision, float valueRangeSize, float medianValue_f)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_4D(oriData, r1, r2, r3, r4, realPrecision);
		updateQuantizationInfo(quantization_intervals);
	}
	else
		quantization_intervals = intvCapacity;

	int i,j,k, reqLength;
	float pred1D, pred2D, pred3D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;

	int dataLength = r1*r2*r3*r4;
	int r234 = r2*r3*r4;
	int r34 = r3*r4;

	P0 = (float*)malloc(r34*sizeof(float));
	P1 = (float*)malloc(r34*sizeof(float));

	float medianValue = medianValue_f;
	short radExpo = getExponent_float(valueRangeSize/2);
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));

	float* spaceFillingValue = oriData; //

	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));

	int l;
	for (l = 0; l < r1; l++)
	{

		///////////////////////////	Process layer-0 ///////////////////////////
		/* Process Row-0 data 0*/
		int index = l*r234;
		int index2D = 0;

		type[index] = 0;
		addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
		compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[index2D] = vce->data;

		/* Process Row-0 data 1*/
		index = l*r234+1;
		index2D = 1;

		pred1D = P1[index2D-1];
		diff = spaceFillingValue[index] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P1[index2D] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
		}
		else
		{
			type[index] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index2D] = vce->data;
		}

		/* Process Row-0 data 2 --> data r4-1 */
		for (j = 2; j < r4; j++)
		{
			index = l*r234+j;
			index2D = j;

			pred1D = 2*P1[index2D-1] - P1[index2D-2];
			diff = spaceFillingValue[index] - pred1D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index2D] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				type[index] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index2D] = vce->data;
			}
		}

		/* Process Row-1 --> Row-r3-1 */
		for (i = 1; i < r3; i++)
		{
			/* Process row-i data 0 */
			index = l*r234+i*r4;
			index2D = i*r4;

			pred1D = P1[index2D-r4];
			diff = spaceFillingValue[index] - pred1D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index2D] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				type[index] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index2D] = vce->data;
			}

			/* Process row-i data 1 --> data r4-1*/
			for (j = 1; j < r4; j++)
			{
				index = l*r234+i*r4+j;
				index2D = i*r4+j;

				pred2D = P1[index2D-1] + P1[index2D-r4] - P1[index2D-r4-1];

				diff = spaceFillingValue[index] - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P1[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				}
				else
				{
					type[index] = 0;

					addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
					compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P1[index2D] = vce->data;
				}
			}
		}


		///////////////////////////	Process layer-1 --> layer-r2-1 ///////////////////////////

		for (k = 1; k < r2; k++)
		{
			/* Process Row-0 data 0*/
			index = l*r234+k*r34;
			index2D = 0;

			pred1D = P1[index2D];
			diff = spaceFillingValue[index] - pred1D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				type[index] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[index2D] = vce->data;
			}

			/* Process Row-0 data 1 --> data r4-1 */
			for (j = 1; j < r4; j++)
			{
				index = l*r234+k*r34+j;
				index2D = j;

				pred2D = P0[index2D-1] + P1[index2D] - P1[index2D-1];
				diff = spaceFillingValue[index] - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				}
				else
				{
					type[index] = 0;

					addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
					compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}
			}

			/* Process Row-1 --> Row-r3-1 */
			for (i = 1; i < r3; i++)
			{
				/* Process Row-i data 0 */
				index = l*r234+k*r34+i*r4;
				index2D = i*r4;

				pred2D = P0[index2D-r4] + P1[index2D] - P1[index2D-r4];
				diff = spaceFillingValue[index] - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
				}
				else
				{
					type[index] = 0;

					addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
					compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}

				/* Process Row-i data 1 --> data r4-1 */
				for (j = 1; j < r4; j++)
				{
					index = l*r234+k*r34+i*r4+j;
					index2D = i*r4+j;

					//printf ("global index = %d, prediction index = %d\n", index, index2D);

					pred3D = P0[index2D-1] + P0[index2D-r4]+ P1[index2D] - P0[index2D-r4-1] - P1[index2D-r4] - P1[index2D-1] + P1[index2D-r4-1];
					diff = spaceFillingValue[index] - pred3D;

					itvNum = fabs(diff)/realPrecision + 1;

					if (itvNum < intvCapacity)
					{
						if (diff < 0) itvNum = -itvNum;
						type[index] = (int) (itvNum/2) + intvRadius;
						P0[index2D] = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
					}
					else
					{
						type[index] = 0;

						addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
						compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
						updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
						memcpy(preDataBytes,vce->curBytes,4);
						addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
						P0[index2D] = vce->data;
					}
				}
			}

			float *Pt;
			Pt = P1;
			P1 = P0;
			P0 = Pt;
		}
	}

	free(P0);
	free(P1);
	int exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitLengthArray->array, resiBitLengthArray->size,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);

	return tdps;
}

void SZ_compress_args_float_NoCkRngeNoGzip_4D(unsigned char** newByteData, float *oriData, int r1, int r2, int r3, int r4, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f)
{
	SZ_Reset();

	TightDataPointStorageF* tdps = SZ_compress_float_4D_MDQ(oriData, r1, r2, r3, r4, realPrecision, valueRangeSize, medianValue_f);

	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);

	int dataLength = r1*r2*r3*r4;
	if(*outSize>dataLength*sizeof(float))
		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);

	free_TightDataPointStorageF(tdps);
}

void SZ_compress_args_float_withinRange(unsigned char** newByteData, float *oriData, int dataLength, int *outSize)
{
	TightDataPointStorageF* tdps = (TightDataPointStorageF*) malloc(sizeof(TightDataPointStorageF));
	tdps->rtypeArray = NULL;
	tdps->typeArray = NULL;	
	tdps->leadNumArray = NULL;
	tdps->residualMidBits = NULL;
	
	tdps->allSameData = 1;
	tdps->dataSeriesLength = dataLength;
	tdps->exactMidBytes = (unsigned char*)malloc(sizeof(unsigned char)*4);
	tdps->pwrErrBoundBytes = NULL;
	tdps->isLossless = 0;
	float value = oriData[0];
	floatToBytes(tdps->exactMidBytes, value);
	tdps->exactMidBytes_size = 4;
	
	int tmpOutSize;
	//unsigned char *tmpByteData;
	convertTDPStoFlatBytes_float(tdps, newByteData, &tmpOutSize);

	//*newByteData = (unsigned char*)malloc(sizeof(unsigned char)*12); //for floating-point data (1+3+4+4)
	//memcpy(*newByteData, tmpByteData, 12);
	*outSize = 12;
	free_TightDataPointStorageF(tdps);	
}

int SZ_compress_args_float_wRngeNoGzip(unsigned char** newByteData, float *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio)
{
	int status = SZ_SCES;
	int dataLength = computeDataLength(r5,r4,r3,r2,r1);
	float valueRangeSize = 0, medianValue = 0;
	
	float min = computeRangeSize_float(oriData, dataLength, &valueRangeSize, &medianValue);
	float max = min+valueRangeSize;
	double realPrecision = getRealPrecision_float(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio, &status);
		
	if(valueRangeSize <= realPrecision)
	{
		SZ_compress_args_float_withinRange(newByteData, oriData, dataLength, outSize);
	}
	else
	{
//		SZ_compress_args_float_NoCkRngeNoGzip_2D(newByteData, oriData, r2, r1, realPrecision, outSize);
		if(r5==0&&r4==0&&r3==0&&r2==0)
		{
			if(errBoundMode==PW_REL)
				SZ_compress_args_float_NoCkRngeNoGzip_1D_pwr(newByteData, oriData, r1, outSize, min, max);
			else
				SZ_compress_args_float_NoCkRngeNoGzip_1D(newByteData, oriData, r1, realPrecision, outSize, valueRangeSize, medianValue);
		}
		else if(r5==0&&r4==0&&r3==0)
		{
			if(errBoundMode==PW_REL)
				SZ_compress_args_float_NoCkRngeNoGzip_2D_pwr(newByteData, oriData, r2, r1, outSize, min, max);
			else
				SZ_compress_args_float_NoCkRngeNoGzip_2D(newByteData, oriData, r2, r1, realPrecision, outSize, valueRangeSize, medianValue);
		}
		else if(r5==0&&r4==0)
			SZ_compress_args_float_NoCkRngeNoGzip_3D(newByteData, oriData, r3, r2, r1, realPrecision, outSize, valueRangeSize, medianValue);
		else if(r5==0)
			SZ_compress_args_float_NoCkRngeNoGzip_3D(newByteData, oriData, r4*r3, r2, r1, realPrecision, outSize, valueRangeSize, medianValue);
	}
	return status;
}

int SZ_compress_args_float(unsigned char** newByteData, float *oriData, 
int r5, int r4, int r3, int r2, int r1, 
int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio)
{
	errorBoundMode = errBoundMode;
	if(errBoundMode==PW_REL)
		pw_relBoundRatio = relBoundRatio;	
	int status = SZ_SCES;
	int dataLength = computeDataLength(r5,r4,r3,r2,r1);
	float valueRangeSize = 0, medianValue = 0;
	
	float min = computeRangeSize_float(oriData, dataLength, &valueRangeSize, &medianValue);
	float max = min+valueRangeSize;
	double realPrecision = getRealPrecision_float(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio, &status);
		
	if(valueRangeSize <= realPrecision)
	{
		SZ_compress_args_float_withinRange(newByteData, oriData, dataLength, outSize);
	}
	else
	{
		int tmpOutSize = 0;
		unsigned char* tmpByteData;
		if (r2==0)
		{
			if(errBoundMode==PW_REL)
				SZ_compress_args_float_NoCkRngeNoGzip_1D_pwr(&tmpByteData, oriData, r1, &tmpOutSize, min, max);
			else
				SZ_compress_args_float_NoCkRngeNoGzip_1D(&tmpByteData, oriData, r1, realPrecision, &tmpOutSize, valueRangeSize, medianValue);
		}
		else
		if (r3==0)
		{
			if(errBoundMode==PW_REL)
				SZ_compress_args_float_NoCkRngeNoGzip_2D_pwr(&tmpByteData, oriData, r2, r1, &tmpOutSize, min, max);
			else
				SZ_compress_args_float_NoCkRngeNoGzip_2D(&tmpByteData, oriData, r2, r1, realPrecision, &tmpOutSize, valueRangeSize, medianValue);
		}
		else
		if (r4==0)
		{
			if(errBoundMode==PW_REL)
				SZ_compress_args_float_NoCkRngeNoGzip_3D_pwr(&tmpByteData, oriData, r3, r2, r1, &tmpOutSize, min, max);
			else
				SZ_compress_args_float_NoCkRngeNoGzip_3D(&tmpByteData, oriData, r3, r2, r1, realPrecision, &tmpOutSize, valueRangeSize, medianValue);
		}
		else
		if (r5==0)
		{
			if(errBoundMode==PW_REL)
				SZ_compress_args_float_NoCkRngeNoGzip_3D_pwr(&tmpByteData, oriData, r4*r3, r2, r1, &tmpOutSize, min, max);
				//ToDO
				//SZ_compress_args_float_NoCkRngeNoGzip_4D_pwr(&tmpByteData, oriData, r4, r3, r2, r1, &tmpOutSize, min, max);
			else
				SZ_compress_args_float_NoCkRngeNoGzip_4D(&tmpByteData, oriData, r4, r3, r2, r1, realPrecision, &tmpOutSize, valueRangeSize, medianValue);
		}
		else
		{
			printf("Error: doesn't support 5 dimensions for now.\n");
			status = SZ_DERR; //dimension error
		}
		//Call Gzip to do the further compression.
		if(szMode==SZ_BEST_SPEED)
		{
			*outSize = tmpOutSize;
			*newByteData = tmpByteData;			
		}
		else if(szMode==SZ_BEST_COMPRESSION || szMode==SZ_DEFAULT_COMPRESSION)
		{
			*outSize = (int)zlib_compress2(tmpByteData, tmpOutSize, newByteData, gzipMode);
			free(tmpByteData);
		}
		else
		{
			printf("Error: Wrong setting of szMode in the float compression.\n");
			status = SZ_MERR; //mode error			
		}
	}
	SZ_ReleaseHuffman();
	return status;
}

/**
 * 
 * 
 * @return status SUCCESSFUL (SZ_SCES) or not (other error codes) 
 * */
int SZ_decompress_args_float(float** newData, int r5, int r4, int r3, int r2, int r1, unsigned char* cmpBytes, int cmpSize)
{
	int status = SZ_SCES;
	int dataLength = computeDataLength(r5,r4,r3,r2,r1);
	
	//unsigned char* tmpBytes;
	int targetUncompressSize = dataLength <<2; //i.e., *4
	//tmpSize must be "much" smaller than dataLength
	int tmpSize = 12, i;
	unsigned char* szTmpBytes;	
	
	if(cmpSize!=12)
	{
		int isZlib = isZlibFormat(cmpBytes[0], cmpBytes[1]);
		if(isZlib)
			szMode = SZ_BEST_COMPRESSION;
		else
			szMode = SZ_BEST_SPEED;		
		if(szMode==SZ_BEST_SPEED)
		{
			tmpSize = cmpSize;
			szTmpBytes = cmpBytes;	
		}
		else if(szMode==SZ_BEST_COMPRESSION || szMode==SZ_DEFAULT_COMPRESSION)
		{
			if(targetUncompressSize<MIN_ZLIB_DEC_ALLOMEM_BYTES) //Considering the minimum size
				targetUncompressSize = MIN_ZLIB_DEC_ALLOMEM_BYTES; 
			tmpSize = zlib_uncompress2(cmpBytes, (unsigned long)cmpSize, &szTmpBytes, (unsigned long)targetUncompressSize);			
			//szTmpBytes = (unsigned char*)malloc(sizeof(unsigned char)*tmpSize);
			//memcpy(szTmpBytes, tmpBytes, tmpSize);
			//free(tmpBytes); //release useless memory		
		}
		else
		{
			printf("Wrong value of szMode in the double compressed bytes.\n");
			status = SZ_MERR;
			return status;
		}	
	}
	else
		szTmpBytes = cmpBytes;
	
	//TODO: convert szTmpBytes to data array.
	TightDataPointStorageF* tdps;
	int errBoundMode = new_TightDataPointStorageF_fromFlatBytes(&tdps, szTmpBytes, tmpSize);
	
	//writeByteData(tdps->typeArray, tdps->typeArray_size, "decompress-typebytes.tbt");
	
	int dim = computeDimension(r5,r4,r3,r2,r1);	
	int floatSize = sizeof(float);
	if(tdps->isLossless)
	{
		*newData = (float*)malloc(floatSize*dataLength);
		if(sysEndianType==BIG_ENDIAN_SYSTEM)
		{
			memcpy(*newData, szTmpBytes+8, dataLength*floatSize);
		}
		else
		{
			unsigned char* p = szTmpBytes+8;
			for(i=0;i<dataLength;i++,p+=floatSize)
				(*newData)[i] = bytesToFloat(p);
		}		
	}
	else if (dim == 1)
		getSnapshotData_float_1D(newData,r1,tdps, errBoundMode);
	else
	if (dim == 2)
		getSnapshotData_float_2D(newData,r2,r1,tdps, errBoundMode);
	else
	if (dim == 3)
		getSnapshotData_float_3D(newData,r3,r2,r1,tdps, errBoundMode);
	else
	if (dim == 4)
		getSnapshotData_float_4D(newData,r4,r3,r2,r1,tdps, errBoundMode);
	else
	{
		printf("Error: currently support only at most 4 dimensions!\n");
		status = SZ_DERR;
	}
	free_TightDataPointStorageF(tdps);
	if(szMode!=SZ_BEST_SPEED && cmpSize!=12)
		free(szTmpBytes);
	SZ_ReleaseHuffman();	
	return status;
}

void computeReqLength_float(float realPrecision, short radExpo, int* reqLength, float* medianValue)
{
	short reqExpo = getPrecisionReqLength_double(realPrecision);
	*reqLength = 9+radExpo - reqExpo; //radExpo-reqExpo == reqMantiLength
	if(*reqLength<9)
		*reqLength = 9;
	if(*reqLength>32)
	{	
		*reqLength = 32;
		*medianValue = 0;
	}			
}

//TODO
int SZ_compress_args_float_subblock(unsigned char* compressedBytes, float *oriData,
int r5, int r4, int r3, int r2, int r1,
int s5, int s4, int s3, int s2, int s1,
int e5, int e4, int e3, int e2, int e1,
int *outSize, int errBoundMode, double absErr_Bound, double relBoundRatio)
{
	int status = SZ_SCES;
	float valueRangeSize = 0, medianValue = 0;
	float min = computeRangeSize_float_subblock(oriData, &valueRangeSize, &medianValue, r5, r4, r3, r2, r1, s5, s4, s3, s2, s1, e5, e4, e3, e2, e1);
	float max = min+valueRangeSize;

	double realPrecision = getRealPrecision_float(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio, &status);

	if(valueRangeSize <= realPrecision)
	{
		//TODO
		//SZ_compress_args_float_withinRange_subblock();
	}
	else
	{
		if (r2==0)
		{
			if(errBoundMode==PW_REL)
			{
				//TODO
				//SZ_compress_args_float_NoCkRngeNoGzip_1D_pwr_subblock();
				printf ("Current subblock version does not support point-wise relative error bound.\n");
			}
			else
				SZ_compress_args_float_NoCkRnge_1D_subblock(compressedBytes, oriData, realPrecision, outSize, valueRangeSize, medianValue, r1, s1, e1);
		}
		else
		if (r3==0)
		{
			//TODO
			if(errBoundMode==PW_REL)
			{
				//TODO
				//SZ_compress_args_float_NoCkRngeNoGzip_2D_pwr_subblock();
				printf ("Current subblock version does not support point-wise relative error bound.\n");
			}
			else
				SZ_compress_args_float_NoCkRnge_2D_subblock(compressedBytes, oriData, realPrecision, outSize, valueRangeSize, medianValue, r2, r1, s2, s1, e2, e1);
		}
		else
		if (r4==0)
		{
			if(errBoundMode==PW_REL)
			{
				//TODO
				//SZ_compress_args_float_NoCkRngeNoGzip_3D_pwr_subblock();
				printf ("Current subblock version does not support point-wise relative error bound.\n");
			}
			else
				SZ_compress_args_float_NoCkRnge_3D_subblock(compressedBytes, oriData, realPrecision, outSize, valueRangeSize, medianValue, r3, r2, r1, s3, s2, s1, e3, e2, e1);
		}
		else
		if (r5==0)
		{
			if(errBoundMode==PW_REL)
			{
				//TODO
				//SZ_compress_args_float_NoCkRngeNoGzip_4D_pwr_subblock();
				printf ("Current subblock version does not support point-wise relative error bound.\n");
			}
			else
				SZ_compress_args_float_NoCkRnge_4D_subblock(compressedBytes, oriData, realPrecision, outSize, valueRangeSize, medianValue, r4, r3, r2, r1, s4, s3, s2, s1, e4, e3, e2, e1);
		}
		else
		{
			printf("Error: doesn't support 5 dimensions for now.\n");
			status = SZ_DERR; //dimension error
		}
	}
	SZ_ReleaseHuffman();
	return status;
}

void SZ_compress_args_float_NoCkRnge_1D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f,
int r1, int s1, int e1)
{
	SZ_Reset();
	TightDataPointStorageF* tdps = SZ_compress_float_1D_MDQ_subblock(oriData, realPrecision, valueRangeSize, medianValue_f, r1, s1, e1);

	if (szMode==SZ_BEST_SPEED)
		convertTDPStoFlatBytes_float_args(tdps, compressedBytes, outSize);
	else if(szMode==SZ_BEST_COMPRESSION || szMode==SZ_DEFAULT_COMPRESSION)
	{
		unsigned char *tmpCompBytes;
		int tmpOutSize;
		convertTDPStoFlatBytes_float(tdps, &tmpCompBytes, &tmpOutSize);
		*outSize = (int)zlib_compress3(tmpCompBytes, tmpOutSize, compressedBytes, gzipMode);
		free(tmpCompBytes);
	}
	else
	{
		printf ("Error: Wrong setting of szMode in the double compression.\n");
	}

	//TODO
//	if(*outSize>dataLength*sizeof(float))
//		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);

	free_TightDataPointStorageF(tdps);
}

void SZ_compress_args_float_NoCkRnge_2D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f,
int r2, int r1, int s2, int s1, int e2, int e1)
{
	SZ_Reset();
	TightDataPointStorageF* tdps = SZ_compress_float_2D_MDQ_subblock(oriData, realPrecision, valueRangeSize, medianValue_f, r2, r1, s2, s1, e2, e1);

	if (szMode==SZ_BEST_SPEED)
		convertTDPStoFlatBytes_float_args(tdps, compressedBytes, outSize);
	else if(szMode==SZ_BEST_COMPRESSION || szMode==SZ_DEFAULT_COMPRESSION)
	{
		unsigned char *tmpCompBytes;
		int tmpOutSize;
		convertTDPStoFlatBytes_float(tdps, &tmpCompBytes, &tmpOutSize);
		*outSize = (int)zlib_compress3(tmpCompBytes, tmpOutSize, compressedBytes, gzipMode);
		free(tmpCompBytes);
	}
	else
	{
		printf ("Error: Wrong setting of szMode in the double compression.\n");
	}

	//TODO
//	if(*outSize>dataLength*sizeof(float))
//		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);

	free_TightDataPointStorageF(tdps);
}

void SZ_compress_args_float_NoCkRnge_3D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f,
int r3, int r2, int r1, int s3, int s2, int s1, int e3, int e2, int e1)
{
	SZ_Reset();
	TightDataPointStorageF* tdps = SZ_compress_float_3D_MDQ_subblock(oriData, realPrecision, valueRangeSize, medianValue_f, r3, r2, r1, s3, s2, s1, e3, e2, e1);

	if (szMode==SZ_BEST_SPEED)
		convertTDPStoFlatBytes_float_args(tdps, compressedBytes, outSize);
	else if(szMode==SZ_BEST_COMPRESSION || szMode==SZ_DEFAULT_COMPRESSION)
	{
		unsigned char *tmpCompBytes;
		int tmpOutSize;
		convertTDPStoFlatBytes_float(tdps, &tmpCompBytes, &tmpOutSize);
		*outSize = (int)zlib_compress3(tmpCompBytes, tmpOutSize, compressedBytes, gzipMode);
		free(tmpCompBytes);
	}
	else
	{
		printf ("Error: Wrong setting of szMode in the double compression.\n");
	}

	//TODO
//	if(*outSize>dataLength*sizeof(float))
//		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);

	free_TightDataPointStorageF(tdps);
}

void SZ_compress_args_float_NoCkRnge_4D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f,
int r4, int r3, int r2, int r1, int s4, int s3, int s2, int s1, int e4, int e3, int e2, int e1)
{
	SZ_Reset();
	TightDataPointStorageF* tdps = SZ_compress_float_4D_MDQ_subblock(oriData, realPrecision, valueRangeSize, medianValue_f, r4, r3, r2, r1, s4, s3, s2, s1, e4, e3, e2, e1);

	if (szMode==SZ_BEST_SPEED)
		convertTDPStoFlatBytes_float_args(tdps, compressedBytes, outSize);
	else if(szMode==SZ_BEST_COMPRESSION || szMode==SZ_DEFAULT_COMPRESSION)
	{
		unsigned char *tmpCompBytes;
		int tmpOutSize;
		convertTDPStoFlatBytes_float(tdps, &tmpCompBytes, &tmpOutSize);
		*outSize = (int)zlib_compress3(tmpCompBytes, tmpOutSize, compressedBytes, gzipMode);
		free(tmpCompBytes);
	}
	else
	{
		printf ("Error: Wrong setting of szMode in the double compression.\n");
	}

	//TODO
//	if(*outSize>dataLength*sizeof(float))
//		SZ_compress_args_float_StoreOriData(oriData, dataLength, tdps, newByteData, outSize);

	free_TightDataPointStorageF(tdps);

}

unsigned int optimize_intervals_float_1D_subblock(float *oriData, double realPrecision, int r1, int s1, int e1)
{
	int dataLength = e1 - s1 + 1;
	oriData = oriData + s1;

	int i = 0;
	unsigned long radiusIndex;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	int totalSampleSize = dataLength/sampleDistance;
	for(i=2;i<dataLength;i++)
	{
		if(i%sampleDistance==0)
		{
			pred_value = 2*oriData[i-1] - oriData[i-2];
			//pred_value = oriData[i-1];
			pred_err = fabs(pred_value - oriData[i]);
			radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
			if(radiusIndex>=maxRangeRadius)
				radiusIndex = maxRangeRadius - 1;
			intervals[radiusIndex]++;
		}
	}
	//compute the appropriate number
	int targetCount = (int)(totalSampleSize*predThreshold);
	int sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;

	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	//printf("accIntervals=%d, powerOf2=%d\n", accIntervals, powerOf2);
	return powerOf2;
}

unsigned int optimize_intervals_float_2D_subblock(float *oriData, double realPrecision, int r1, int r2, int s1, int s2, int e1, int e2)
{
	int R1 = e1 - s1 + 1;
	int R2 = e2 - s2 + 1;

	int i,j, index;
	unsigned long radiusIndex;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	int totalSampleSize = R1*R2/sampleDistance;
	for(i=s1;i<=e1;i++)
	{
		for(j=s2;j<=e2;j++)
		{
			if((i+j)%sampleDistance==0)
			{
				index = i*r2+j;
				pred_value = oriData[index-1] + oriData[index-r2] - oriData[index-r2-1];
				pred_err = fabs(pred_value - oriData[index]);
				radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
				if(radiusIndex>=maxRangeRadius)
					radiusIndex = maxRangeRadius - 1;
				intervals[radiusIndex]++;
			}
		}
	}
	//compute the appropriate number
	int targetCount = (int)(totalSampleSize*predThreshold);
	int sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	//printf("maxRangeRadius = %d, accIntervals=%d, powerOf2=%d\n", maxRangeRadius, accIntervals, powerOf2);
	return powerOf2;
}

unsigned int optimize_intervals_float_3D_subblock(float *oriData, double realPrecision, int r1, int r2, int r3, int s1, int s2, int s3, int e1, int e2, int e3)
{
	int R1 = e1 - s1 + 1;
	int R2 = e2 - s2 + 1;
	int R3 = e3 - s3 + 1;

	int r23 = r2*r3;

	int i,j,k, index;
	unsigned long radiusIndex;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	int totalSampleSize = R1*R2*R3/sampleDistance;
	for(i=s1;i<=e1;i++)
	{
		for(j=s2;j<=e2;j++)
		{
			for(k=s3;k<=e3;k++)
			{
				if((i+j+k)%sampleDistance==0)
				{
					index = i*r23+j*r3+k;
					pred_value = oriData[index-1] + oriData[index-r3] + oriData[index-r23]
					- oriData[index-1-r23] - oriData[index-r3-1] - oriData[index-r3-r23] + oriData[index-r3-r23-1];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
					if(radiusIndex>=maxRangeRadius)
						radiusIndex = maxRangeRadius - 1;
					intervals[radiusIndex]++;
				}
			}
		}
	}
	//compute the appropriate number
	int targetCount = (int)(totalSampleSize*predThreshold);
	int sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	return powerOf2;
}

unsigned int optimize_intervals_float_4D_subblock(float *oriData, double realPrecision,
int r1, int r2, int r3, int r4, int s1, int s2, int s3, int s4, int e1, int e2, int e3, int e4)
{
	int R1 = e1 - s1 + 1;
	int R2 = e2 - s2 + 1;
	int R3 = e3 - s3 + 1;
	int R4 = e4 - s4 + 1;

	int r34 = r3*r4;
	int r234 = r2*r3*r4;

	int i,j,k,l, index;
	unsigned long radiusIndex;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	int totalSampleSize = R1*R2*R3*R4/sampleDistance;
	for(i=s1;i<=e1;i++)
	{
		for(j=s2;j<=e2;j++)
		{
			for(k=s3;k<=e3;k++)
			{
				for (l=s4;l<=e4;l++)
				{
					if((i+j+k+l)%sampleDistance==0)
					{
						index = i*r234+j*r34+k*r4;
						pred_value = oriData[index-1] + oriData[index-r4] + oriData[index-r34]
									- oriData[index-1-r34] - oriData[index-r4-1] - oriData[index-r4-r34] + oriData[index-r4-r34-1];
						pred_err = fabs(pred_value - oriData[index]);
						radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
						if(radiusIndex>=maxRangeRadius)
							radiusIndex = maxRangeRadius - 1;
						intervals[radiusIndex]++;
					}
				}
			}
		}
	}
	//compute the appropriate number
	int targetCount = (int)(totalSampleSize*predThreshold);
	int sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	return powerOf2;
}

TightDataPointStorageF* SZ_compress_float_1D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
int r1, int s1, int e1)
{
	int dataLength = e1 - s1 + 1;
	unsigned int quantization_intervals;
	if(optQuantMode==1)
		quantization_intervals = optimize_intervals_float_1D_subblock(oriData, realPrecision, r1, s1, e1);
	else
		quantization_intervals = intvCapacity;
	updateQuantizationInfo(quantization_intervals);

	int i, reqLength;
	float medianValue = medianValue_f;
	short reqExpo = getPrecisionReqLength_float((float)realPrecision);
	short radExpo = getExponent_float(valueRangeSize/2);

	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));

	float* spaceFillingValue = oriData + s1;

	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	type[0] = 0;

	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;
	float last3CmprsData[3] = {0};

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));

	//add the first data
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_float(last3CmprsData, vce->data);

	//add the second data
	type[1] = 0;
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
	compressSingleFloatValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_float(last3CmprsData, vce->data);

	int state;
	float lcf, qcf;
	double checkRadius;
	float curData;
	float pred;
	float predAbsErr;
	float min_pred, minErr, minIndex;
	int a = 0;
	checkRadius = (intvCapacity-1)*realPrecision;
	double interval = 2*realPrecision;

	for(i=2;i<dataLength;i++)
	{
		curData = spaceFillingValue[i];
		pred = 2*last3CmprsData[0] - last3CmprsData[1];
		predAbsErr = fabs(curData - pred);
		if(predAbsErr<=checkRadius)
		{
			state = (predAbsErr/realPrecision+1)/2;
			if(curData>=pred)
			{
				type[i] = intvRadius+state;
				pred = pred + state*interval;
			}
			else
			{
				type[i] = intvRadius-state;
				pred = pred - state*interval;
			}

			listAdd_float(last3CmprsData, pred);
			continue;
		}

		//unpredictable data processing
		type[i] = 0;
		addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);

		compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);

		listAdd_float(last3CmprsData, vce->data);
	}

	int exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitLengthArray->array, resiBitLengthArray->size,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);

	return tdps;
}

TightDataPointStorageF* SZ_compress_float_2D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
int r1, int r2, int s1, int s2, int e1, int e2)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_2D_subblock(oriData, realPrecision, r1, r2, s1, s2, e1, e2);
		updateQuantizationInfo(quantization_intervals);
	}
	else
		quantization_intervals = intvCapacity;

	int i,j, reqLength;
	float pred1D, pred2D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;

	int R1 = e1 - s1 + 1;
	int R2 = e2 - s2 + 1;
	int dataLength = R1*R2;

	P0 = (float*)malloc(R2*sizeof(float));
	memset(P0, 0, R2*sizeof(float));
	P1 = (float*)malloc(R2*sizeof(float));
	memset(P1, 0, R2*sizeof(float));

	float medianValue = medianValue_f;
	short radExpo = getExponent_float(valueRangeSize/2);
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));

	float* spaceFillingValue = oriData; //

	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));

	/* Process Row-s1 data s2*/
	int gIndex;
	int lIndex;

	gIndex = s1*r2+s2;
	lIndex = 0;

	type[lIndex] = 0;
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
	compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[0] = vce->data;

	/* Process Row-s1 data s2+1*/
	gIndex = s1*r2+(s2+1);
	lIndex = 1;

	pred1D = P1[0];
	diff = spaceFillingValue[gIndex] - pred1D;

	itvNum =  fabs(diff)/realPrecision + 1;

	if (itvNum < intvCapacity)
	{
		if (diff < 0) itvNum = -itvNum;
		type[lIndex] = (int) (itvNum/2) + intvRadius;
		P1[1] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
	}
	else
	{
		type[lIndex] = 0;

		addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
		compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[1] = vce->data;
	}

    /* Process Row-s1 data s2+2 --> data e2 */
	for (j = 2; j < R2; j++)
	{
		gIndex = s1*r2+(s2+j);
		lIndex = j;

		pred1D = 2*P1[j-1] - P1[j-2];
		diff = spaceFillingValue[gIndex] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[lIndex] = (int) (itvNum/2) + intvRadius;
			P1[j] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
		}
		else
		{
			type[lIndex] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[j] = vce->data;
		}
	}

	/* Process Row-s1+1 --> Row-e1 */
	for (i = 1; i < R1; i++)
	{
		/* Process row-s1+i data s2 */
		gIndex = (s1+i)*r2+s2;
		lIndex = i*R2;

		pred1D = P1[0];
		diff = spaceFillingValue[gIndex] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[lIndex] = (int) (itvNum/2) + intvRadius;
			P0[0] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
		}
		else
		{
			type[lIndex] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[0] = vce->data;
		}

		/* Process row-s1+i data s2+1 --> e2 */
		for (j = 1; j < R2; j++)
		{
			gIndex = (s1+i)*r2+(s2+j);
			lIndex = i*R2+j;

//			printf ("global index = %d, local index = %d\n", gIndex, lIndex);

			pred2D = P0[j-1] + P1[j] - P1[j-1];

			diff = spaceFillingValue[gIndex] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P0[j] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[j] = vce->data;
			}
		}

		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}

	free(P0);
	free(P1);
	int exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitLengthArray->array, resiBitLengthArray->size,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);

	return tdps;
}

TightDataPointStorageF* SZ_compress_float_3D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
int r1, int r2, int r3, int s1, int s2, int s3, int e1, int e2, int e3)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_3D_subblock(oriData, realPrecision, r1, r2, r3, s1, s2, s3, e1, e2, e3);
		updateQuantizationInfo(quantization_intervals);
	}
	else
		quantization_intervals = intvCapacity;

	int i,j,k, reqLength;
	float pred1D, pred2D, pred3D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;

	int R1 = e1 - s1 + 1;
	int R2 = e2 - s2 + 1;
	int R3 = e3 - s3 + 1;
	int dataLength = R1*R2*R3;

	int r23 = r2*r3;
	int R23 = R2*R3;

	P0 = (float*)malloc(R23*sizeof(float));
	P1 = (float*)malloc(R23*sizeof(float));

	float medianValue = medianValue_f;
	short radExpo = getExponent_float(valueRangeSize/2);
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));
	//type[dataLength]=0;

	float* spaceFillingValue = oriData; //

	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));


	///////////////////////////	Process layer-s1 ///////////////////////////
	/* Process Row-s2 data s3*/
	int gIndex; 	//global index
	int lIndex; 	//local index
	int index2D; 	//local 2D index

	gIndex = s1*r23+s2*r3+s3;
	lIndex = 0;
	index2D = 0;

	type[lIndex] = 0;
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
	compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[index2D] = vce->data;

	/* Process Row-s2 data s3+1*/
	gIndex = s1*r23+s2*r3+s3+1;
	lIndex = 1;
	index2D = 1;

	pred1D = P1[index2D-1];
	diff = spaceFillingValue[gIndex] - pred1D;

	itvNum = fabs(diff)/realPrecision + 1;

	if (itvNum < intvCapacity)
	{
		if (diff < 0) itvNum = -itvNum;
		type[lIndex] = (int) (itvNum/2) + intvRadius;
		P1[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
	}
	else
	{
		type[lIndex] = 0;

		addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
		compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[index2D] = vce->data;
	}

    /* Process Row-s2 data s3+2 --> data e3 */
	for (j = 2; j < R3; j++)
	{
		gIndex = s1*r23+s2*r3+s3+j;
		lIndex = j;
		index2D = j;

		pred1D = 2*P1[index2D-1] - P1[index2D-2];
		diff = spaceFillingValue[gIndex] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[lIndex] = (int) (itvNum/2) + intvRadius;
			P1[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
		}
		else
		{
			type[lIndex] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index2D] = vce->data;
		}
	}

	/* Process Row-s2+1 --> Row-e2 */
	for (i = 1; i < R2; i++)
	{
		/* Process row-s2+i data s3 */
		gIndex = s1*r23+(s2+i)*r3+s3;
		lIndex = i*R3;
		index2D = i*R3;

		pred1D  = P1[index2D-R3];
		diff = spaceFillingValue[gIndex] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[lIndex] = (int) (itvNum/2) + intvRadius;
			P1[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
		}
		else
		{
			type[lIndex] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index2D] = vce->data;
		}

		/* Process row-s2+i data s3+1 --> data e3*/
		for (j = 1; j < R3; j++)
		{
			gIndex = s1*r23+(s2+i)*r3+s3+j;
			lIndex = i*R3+j;
			index2D = i*R3+j;

			pred2D  = P1[index2D-1] + P1[index2D-R3] - P1[index2D-R3-1];
			diff = spaceFillingValue[gIndex] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P1[index2D] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index2D] = vce->data;
			}
		}
	}


	///////////////////////////	Process layer-s1+1 --> layer-e1 ///////////////////////////

	for (k = 1; k < R1; k++)
	{
		/* Process Row-s2 data s3*/
		gIndex = (s1+k)*r23+s2*r3+s3;
		lIndex = k*R23;
		index2D = 0;

		pred1D = P1[index2D];
		diff = spaceFillingValue[gIndex] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[lIndex] = (int) (itvNum/2) + intvRadius;
			P0[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
		}
		else
		{
			type[lIndex] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[index2D] = vce->data;
		}

	    /* Process Row-s2 data s3+1 --> data e3 */
		for (j = 1; j < R3; j++)
		{
			gIndex = (s1+k)*r23+s2*r3+s3+j;
			lIndex = k*R23+j;
			index2D = j;

			pred2D = P0[index2D-1] + P1[index2D] - P1[index2D-1];
			diff = spaceFillingValue[gIndex] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[index2D] = vce->data;
			}
		}

	    /* Process Row-s2+1 --> Row-e2 */
		for (i = 1; i < R2; i++)
		{
			/* Process Row-s2+i data s3 */
			gIndex = (s1+k)*r23+(s2+i)*r3+s3;
			lIndex = k*R23+i*R3;
			index2D = i*R3;

			pred2D = P0[index2D-R3] + P1[index2D] - P1[index2D-R3];
			diff = spaceFillingValue[gIndex] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[index2D] = vce->data;
			}

			/* Process Row-s2+i data s3+1 --> data e3 */
			for (j = 1; j < R3; j++)
			{
				gIndex = (s1+k)*r23+(s2+i)*r3+s3+j;
				lIndex = k*R23+i*R3+j;
				index2D = i*R3+j;

//				printf ("global index = %d, local index = %d\n", gIndex, lIndex);

				pred3D = P0[index2D-1] + P0[index2D-R3]+ P1[index2D] - P0[index2D-R3-1] - P1[index2D-R3] - P1[index2D-1] + P1[index2D-R3-1];
				diff = spaceFillingValue[gIndex] - pred3D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[lIndex] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred3D + 2 * (type[lIndex] - intvRadius) * realPrecision;
				}
				else
				{
					type[lIndex] = 0;

					addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
					compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}
			}
		}

		float *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}

	free(P0);
	free(P1);
	int exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitLengthArray->array, resiBitLengthArray->size,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);

	return tdps;
}

TightDataPointStorageF* SZ_compress_float_4D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
int r1, int r2, int r3, int r4, int s1, int s2, int s3, int s4, int e1, int e2, int e3, int e4)
{
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_4D_subblock(oriData, realPrecision, r1, r2, r3, r4, s1, s2, s3, s4, e1, e2, e3, e4);
		updateQuantizationInfo(quantization_intervals);
	}
	else
		quantization_intervals = intvCapacity;

	int i,j,k, reqLength;
	float pred1D, pred2D, pred3D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;

	int R1 = e1 - s1 + 1;
	int R2 = e2 - s2 + 1;
	int R3 = e3 - s3 + 1;
	int R4 = e4 - s4 + 1;

	int dataLength = R1*R2*R3*R4;

	int r34 = r3*r4;
	int r234 = r2*r3*r4;
	int R34 = R3*R4;
	int R234 = R2*R3*R4;

	P0 = (float*)malloc(R34*sizeof(float));
	P1 = (float*)malloc(R34*sizeof(float));

	float medianValue = medianValue_f;
	short radExpo = getExponent_float(valueRangeSize/2);
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));

	float* spaceFillingValue = oriData; //

	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	unsigned char preDataBytes[4];
	intToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));


	int l;
	for (l = 0; l < R1; l++)
	{

		///////////////////////////	Process layer-s2 ///////////////////////////
		/* Process Row-s3 data s4*/
		int gIndex; 	//global index
		int lIndex; 	//local index
		int index2D; 	//local 2D index

		gIndex = (s1+l)*r234+s2*r34+s3*r4+s4;
		lIndex = l*R234;
		index2D = 0;

		type[lIndex] = 0;
		addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
		compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[index2D] = vce->data;

		/* Process Row-s3 data s4+1*/
		gIndex = (s1+l)*r234+s2*r34+s3*r4+s4+1;
		lIndex = l*R234+1;
		index2D = 1;

		pred1D = P1[index2D-1];
		diff = spaceFillingValue[gIndex] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[lIndex] = (int) (itvNum/2) + intvRadius;
			P1[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
		}
		else
		{
			type[lIndex] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index2D] = vce->data;
		}

		/* Process Row-s3 data s4+2 --> data e4 */
		for (j = 2; j < R4; j++)
		{
			gIndex = (s1+l)*r234+s2*r34+s3*r4+s4+j;
			lIndex = l*R234+j;
			index2D = j;

			pred1D = 2*P1[index2D-1] - P1[index2D-2];
			diff = spaceFillingValue[gIndex] - pred1D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P1[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index2D] = vce->data;
			}
		}

		/* Process Row-s3+1 --> Row-e3 */
		for (i = 1; i < R3; i++)
		{
			/* Process row-s2+i data s3 */
			gIndex = (s1+l)*r234+s2*r34+(s3+i)*r4+s4;
			lIndex = l*R234+i*R4;
			index2D = i*R4;

			pred1D  = P1[index2D-R4];
			diff = spaceFillingValue[gIndex] - pred1D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P1[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index2D] = vce->data;
			}

			/* Process row-s3+i data s4+1 --> data e4*/
			for (j = 1; j < R4; j++)
			{
				gIndex = (s1+l)*r234+s2*r34+(s3+i)*r4+s4+j;
				lIndex = l*R234+i*R4+j;
				index2D = i*R4+j;

				pred2D  = P1[index2D-1] + P1[index2D-R4] - P1[index2D-R4-1];
				diff = spaceFillingValue[gIndex] - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[lIndex] = (int) (itvNum/2) + intvRadius;
					P1[index2D] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
				}
				else
				{
					type[lIndex] = 0;

					addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
					compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P1[index2D] = vce->data;
				}
			}
		}


		///////////////////////////	Process layer-s2+1 --> layer-e2 ///////////////////////////

		for (k = 1; k < R2; k++)
		{
			/* Process Row-s3 data s4*/
			gIndex = (s1+l)*r234+(s2+k)*r34+s3*r4+s4;
			lIndex = l*R234+k*R34;
			index2D = 0;

			pred1D = P1[index2D];
			diff = spaceFillingValue[gIndex] - pred1D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[lIndex] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred1D + 2 * (type[lIndex] - intvRadius) * realPrecision;
			}
			else
			{
				type[lIndex] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[index2D] = vce->data;
			}

			/* Process Row-s3 data s4+1 --> data e4 */
			for (j = 1; j < R4; j++)
			{
				gIndex = (s1+l)*r234+(s2+k)*r34+s3*r4+s4+j;
				lIndex = l*R234+k*R34+j;
				index2D = j;

				pred2D = P0[index2D-1] + P1[index2D] - P1[index2D-1];
				diff = spaceFillingValue[gIndex] - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[lIndex] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
				}
				else
				{
					type[lIndex] = 0;

					addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
					compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}
			}

			/* Process Row-s3+1 --> Row-e3 */
			for (i = 1; i < R3; i++)
			{
				/* Process Row-s3+i data s4 */
				gIndex = (s1+l)*r234+(s2+k)*r34+(s3+i)*r4+s4;
				lIndex = l*R234+k*R34+i*R4;
				index2D = i*R4;

				pred2D = P0[index2D-R4] + P1[index2D] - P1[index2D-R4];
				diff = spaceFillingValue[gIndex] - pred2D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[lIndex] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred2D + 2 * (type[lIndex] - intvRadius) * realPrecision;
				}
				else
				{
					type[lIndex] = 0;

					addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
					compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,4);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}

				/* Process Row-s3+i data s4+1 --> data e4 */
				for (j = 1; j < R4; j++)
				{
					gIndex = (s1+l)*r234+(s2+k)*r34+(s3+i)*r4+s4+j;
					lIndex = l*R234+k*R34+i*R4+j;
					index2D = i*R4+j;

//					printf ("global index = %d, local index = %d\n", gIndex, lIndex);

					pred3D = P0[index2D-1] + P0[index2D-R4]+ P1[index2D] - P0[index2D-R4-1] - P1[index2D-R4] - P1[index2D-1] + P1[index2D-R4-1];
					diff = spaceFillingValue[gIndex] - pred3D;

					itvNum = fabs(diff)/realPrecision + 1;

					if (itvNum < intvCapacity)
					{
						if (diff < 0) itvNum = -itvNum;
						type[lIndex] = (int) (itvNum/2) + intvRadius;
						P0[index2D] = pred3D + 2 * (type[lIndex] - intvRadius) * realPrecision;
					}
					else
					{
						type[lIndex] = 0;

						addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
						compressSingleFloatValue(vce, spaceFillingValue[gIndex], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
						updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
						memcpy(preDataBytes,vce->curBytes,4);
						addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
						P0[index2D] = vce->data;
					}
				}
			}

			float *Pt;
			Pt = P1;
			P1 = P0;
			P0 = Pt;
		}

	}

	free(P0);
	free(P1);
	int exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageF* tdps;

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitLengthArray->array, resiBitLengthArray->size,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	free(vce);
	free(lce);
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);

	return tdps;
}

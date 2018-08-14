/**
 *  @file sz_float_pwr.c
 *  @author Sheng Di
 *  @date Aug, 2016
 *  @brief SZ_Init, Compression and Decompression functions
 * This file contains the compression/decompression functions related to point-wise relative errors
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "sz.h"
//#include "ExpSegmentConstructor.h"
#include "CompressElement.h"
#include "DynamicByteArray.h"
#include "DynamicIntArray.h"
#include "TightDataPointStorageF.h"
#include "sz_float.h"
#include "sz_float_pwr.h"
#include "zlib.h"
#include "rw.h"

void compute_segment_precisions_float_1D(float *oriData, size_t dataLength, float* pwrErrBound, unsigned char* pwrErrBoundBytes, double globalPrecision)
{
	size_t i = 0, j = 0, k = 0;
	float realPrecision = oriData[0]!=0?fabs(pw_relBoundRatio*oriData[0]):pw_relBoundRatio; 
	float approxPrecision;
	unsigned char realPrecBytes[4];
	float curPrecision;
	float curValue;
	float sum = 0;
	for(i=0;i<dataLength;i++)
	{
		curValue = oriData[i];
		if(i%segment_size==0&&i>0)
		{
			//get two first bytes of the realPrecision
			if(pwr_type==SZ_PWR_AVG_TYPE)
			{
				realPrecision = sum/segment_size;
				sum = 0;			
			}
			if(errorBoundMode==ABS_AND_PW_REL||errorBoundMode==REL_AND_PW_REL)
				realPrecision = realPrecision<globalPrecision?realPrecision:globalPrecision; 
			else if(errorBoundMode==ABS_OR_PW_REL||errorBoundMode==REL_OR_PW_REL)
				realPrecision = realPrecision<globalPrecision?globalPrecision:realPrecision;
				
			floatToBytes(realPrecBytes, realPrecision);
			realPrecBytes[2] = realPrecBytes[3] = 0;
			approxPrecision = bytesToFloat(realPrecBytes);
			//put the realPrecision in float* pwrErBound
			pwrErrBound[j++] = approxPrecision;
			//put the two bytes in pwrErrBoundBytes
			pwrErrBoundBytes[k++] = realPrecBytes[0];
			pwrErrBoundBytes[k++] = realPrecBytes[1];
			
			realPrecision = fabs(pw_relBoundRatio*curValue);
		}
		
		if(curValue!=0)
		{
			curPrecision = fabs(pw_relBoundRatio*curValue);
			
			switch(pwr_type)
			{
			case SZ_PWR_MIN_TYPE: 
				if(realPrecision>curPrecision)
					realPrecision = curPrecision;	
				break;
			case SZ_PWR_AVG_TYPE:
				sum += curPrecision;
				break;
			case SZ_PWR_MAX_TYPE:
				if(realPrecision<curPrecision)
					realPrecision = curPrecision;					
				break;
			}
		}
	}
	if(pwr_type==SZ_PWR_AVG_TYPE)
	{
		int size = dataLength%segment_size==0?segment_size:dataLength%segment_size;
		realPrecision = sum/size;		
	}	
	if(errorBoundMode==ABS_AND_PW_REL||errorBoundMode==REL_AND_PW_REL)
		realPrecision = realPrecision<globalPrecision?realPrecision:globalPrecision; 
	else if(errorBoundMode==ABS_OR_PW_REL||errorBoundMode==REL_OR_PW_REL)
		realPrecision = realPrecision<globalPrecision?globalPrecision:realPrecision;	
	floatToBytes(realPrecBytes, realPrecision);
	realPrecBytes[2] = realPrecBytes[3] = 0;
	approxPrecision = bytesToFloat(realPrecBytes);
	//put the realPrecision in float* pwrErBound
	pwrErrBound[j++] = approxPrecision;
	//put the two bytes in pwrErrBoundBytes
	pwrErrBoundBytes[k++] = realPrecBytes[0];
	pwrErrBoundBytes[k++] = realPrecBytes[1];
}

unsigned int optimize_intervals_float_1D_pwr(float *oriData, size_t dataLength, float* pwrErrBound)
{	
	size_t i = 0, j = 0;
	float realPrecision = pwrErrBound[j++];	
	unsigned long radiusIndex;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	int totalSampleSize = dataLength/sampleDistance;
	for(i=2;i<dataLength;i++)
	{
		if(i%segment_size==0)
			realPrecision = pwrErrBound[j++];
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
	size_t targetCount = totalSampleSize*predThreshold;
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

void compute_segment_precisions_float_2D(float *oriData, float* pwrErrBound, 
size_t r1, size_t r2, size_t R2, size_t edgeSize, unsigned char* pwrErrBoundBytes, float Min, float Max, double globalPrecision)
{
	size_t i = 0, j = 0, k = 0, p = 0, index = 0, J; //I=-1,J=-1 if they are needed
	float realPrecision; 
	float approxPrecision;
	unsigned char realPrecBytes[4];
	float curValue, curAbsValue;
	float* statAbsValues = (float*)malloc(R2*sizeof(float));
	
	float max = fabs(Min)<fabs(Max)?fabs(Max):fabs(Min); //get the max abs value.
	float min = fabs(Min)<fabs(Max)?fabs(Min):fabs(Max);
	for(i=0;i<R2;i++)
	{
		if(pwr_type == SZ_PWR_MIN_TYPE)
			statAbsValues[i] = max;
		else if(pwr_type == SZ_PWR_MAX_TYPE)
			statAbsValues[i] = min;
		else
			statAbsValues[i] = 0; //for SZ_PWR_AVG_TYPE
	}
	for(i=0;i<r1;i++)
	{
		for(j=0;j<r2;j++)
		{
			index = i*r2+j;
			curValue = oriData[index];				
			if(((i%edgeSize==edgeSize-1 || i==r1-1) &&j%edgeSize==0&&j>0) || (i%edgeSize==0&&j==0&&i>0))
			{
				if(pwr_type==SZ_PWR_AVG_TYPE)
				{
					int a = edgeSize, b = edgeSize;
					if(j==0)
					{
						if(r2%edgeSize==0) 
							b = edgeSize;
						else
							b = r2%edgeSize;
					}
					if(i==r1-1)
					{
						if(r1%edgeSize==0)
							a = edgeSize;
						else
							a = r1%edgeSize;
					}
					realPrecision = pw_relBoundRatio*statAbsValues[J]/(a*b);
				}
				else
					realPrecision = pw_relBoundRatio*statAbsValues[J];

				if(errorBoundMode==ABS_AND_PW_REL||errorBoundMode==REL_AND_PW_REL)
					realPrecision = realPrecision<globalPrecision?realPrecision:globalPrecision; 
				else if(errorBoundMode==ABS_OR_PW_REL||errorBoundMode==REL_OR_PW_REL)
					realPrecision = realPrecision<globalPrecision?globalPrecision:realPrecision;
					
				floatToBytes(realPrecBytes, realPrecision);
				realPrecBytes[2] = realPrecBytes[3] = 0;
				approxPrecision = bytesToFloat(realPrecBytes);
				//put the realPrecision in float* pwrErBound		
				pwrErrBound[p++] = approxPrecision;
				//put the two bytes in pwrErrBoundBytes
				pwrErrBoundBytes[k++] = realPrecBytes[0];
				pwrErrBoundBytes[k++] = realPrecBytes[1];	
				
				if(pwr_type == SZ_PWR_MIN_TYPE)
					statAbsValues[J] = max;
				else if(pwr_type == SZ_PWR_MAX_TYPE)
					statAbsValues[J] = min;
				else
					statAbsValues[J] = 0; //for SZ_PWR_AVG_TYPE	
			}	
			if(j==0)
				J = 0;
			else if(j%edgeSize==0)
				J++;			
			if(curValue!=0)
			{
				curAbsValue = fabs(curValue);
				
				switch(pwr_type)
				{
				case SZ_PWR_MIN_TYPE: 
					if(statAbsValues[J]>curAbsValue)
						statAbsValues[J] = curAbsValue;	
					break;
				case SZ_PWR_AVG_TYPE:
					statAbsValues[J] += curAbsValue;
					break;
				case SZ_PWR_MAX_TYPE:
					if(statAbsValues[J]<curAbsValue)
						statAbsValues[J] = curAbsValue;					
					break;
				}
			}
		}
	}
		
	if(pwr_type==SZ_PWR_AVG_TYPE)
	{
		int a = edgeSize, b = edgeSize;
		if(r2%edgeSize==0) 
			b = edgeSize;
		else
			b = r2%edgeSize;
		if(r1%edgeSize==0)
			a = edgeSize;
		else
			a = r1%edgeSize;
		realPrecision = pw_relBoundRatio*statAbsValues[J]/(a*b);
	}
	else
		realPrecision = pw_relBoundRatio*statAbsValues[J];		

	if(errorBoundMode==ABS_AND_PW_REL||errorBoundMode==REL_AND_PW_REL)
		realPrecision = realPrecision<globalPrecision?realPrecision:globalPrecision; 
	else if(errorBoundMode==ABS_OR_PW_REL||errorBoundMode==REL_OR_PW_REL)
		realPrecision = realPrecision<globalPrecision?globalPrecision:realPrecision;
		
	floatToBytes(realPrecBytes, realPrecision);
	realPrecBytes[2] = realPrecBytes[3] = 0;
	approxPrecision = bytesToFloat(realPrecBytes);
	//put the realPrecision in float* pwrErBound
	pwrErrBound[p++] = approxPrecision;
	//put the two bytes in pwrErrBoundBytes
	pwrErrBoundBytes[k++] = realPrecBytes[0];
	pwrErrBoundBytes[k++] = realPrecBytes[1];	
	
	free(statAbsValues);
}

unsigned int optimize_intervals_float_2D_pwr(float *oriData, size_t r1, size_t r2, size_t R2, size_t edgeSize, float* pwrErrBound)
{	
	size_t i = 0,j = 0, index, I=0, J=0;
	float realPrecision = pwrErrBound[0];	
	unsigned long radiusIndex;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	size_t totalSampleSize = (r1-1)*(r2-1)/sampleDistance;
	size_t ir2;
	for(i=1;i<r1;i++)
	{
		ir2 = i*r2;
		if(i%edgeSize==0)
		{	
			I++;
			J = 0;
		}
		for(j=1;j<r2;j++)
		{
			index = ir2+j;
			if(j%edgeSize==0)
				J++;
				
			if((i+j)%sampleDistance==0)
			{
				realPrecision = pwrErrBound[I*R2+J];
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
	size_t targetCount = totalSampleSize*predThreshold;
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

void compute_segment_precisions_float_3D(float *oriData, float* pwrErrBound, 
size_t r1, size_t r2, size_t r3, size_t R2, size_t R3, size_t edgeSize, unsigned char* pwrErrBoundBytes, float Min, float Max, double globalPrecision)
{
	size_t i = 0, j = 0, k = 0, p = 0, q = 0, index = 0, J = 0, K = 0; //I=-1,J=-1 if they are needed
	size_t r23 = r2*r3, ir, jr;
	float realPrecision; 
	float approxPrecision;
	unsigned char realPrecBytes[4];
	float curValue, curAbsValue;
	
	float** minAbsValues = create2DArray_float(R2, R3);
	float max = fabs(Min)<fabs(Max)?fabs(Max):fabs(Min); //get the max abs value.	
	for(i=0;i<R2;i++)
		for(j=0;j<R3;j++)
			minAbsValues[i][j] = max;			
			
	for(i=0;i<r1;i++)
	{
		ir = i*r23;		
		if(i%edgeSize==0&&i>0)
		{
			realPrecision = pw_relBoundRatio*minAbsValues[J][K];
			floatToBytes(realPrecBytes, realPrecision);
			memset(&realPrecBytes[2], 0, 2);
			approxPrecision = bytesToFloat(realPrecBytes);
			//put the realPrecision in float* pwrErBound
			pwrErrBound[p++] = approxPrecision;
			//put the two bytes in pwrErrBoundBytes
			//printf("q=%d, i=%d, j=%d, k=%d\n",q,i,j,k);
			pwrErrBoundBytes[q++] = realPrecBytes[0];
			pwrErrBoundBytes[q++] = realPrecBytes[1];
			minAbsValues[J][K] = max;			
		}		
		for(j=0;j<r2;j++)
		{
			jr = j*r3;
			if((i%edgeSize==edgeSize-1 || i == r1-1)&&j%edgeSize==0&&j>0)
			{
				realPrecision = pw_relBoundRatio*minAbsValues[J][K];
				floatToBytes(realPrecBytes, realPrecision);
				memset(&realPrecBytes[2], 0, 2);
				approxPrecision = bytesToFloat(realPrecBytes);
				//put the realPrecision in float* pwrErBound
				pwrErrBound[p++] = approxPrecision;
				//put the two bytes in pwrErrBoundBytes
				//printf("q=%d, i=%d, j=%d, k=%d\n",q,i,j,k);
				pwrErrBoundBytes[q++] = realPrecBytes[0];
				pwrErrBoundBytes[q++] = realPrecBytes[1];
				minAbsValues[J][K] = max;				
			}
			
			if(j==0)
				J = 0;
			else if(j%edgeSize==0)
				J++;					
			
			for(k=0;k<r3;k++)
			{
				index = ir+jr+k;				
				curValue = oriData[index];				
				if((i%edgeSize==edgeSize-1 || i == r1-1)&&(j%edgeSize==edgeSize-1||j==r2-1)&&k%edgeSize==0&&k>0)
				{
					realPrecision = pw_relBoundRatio*minAbsValues[J][K];
					floatToBytes(realPrecBytes, realPrecision);
					memset(&realPrecBytes[2], 0, 2);
					approxPrecision = bytesToFloat(realPrecBytes);
					//put the realPrecision in float* pwrErBound
					pwrErrBound[p++] = approxPrecision;
					//put the two bytes in pwrErrBoundBytes
					//printf("q=%d, i=%d, j=%d, k=%d\n",q,i,j,k);
					pwrErrBoundBytes[q++] = realPrecBytes[0];
					pwrErrBoundBytes[q++] = realPrecBytes[1];
					minAbsValues[J][K] = max;
				}	

				if(k==0)
					K = 0;
				else if(k%edgeSize==0)
					K++;
					
				if(curValue!=0)
				{
					curAbsValue = fabs(curValue);
					if(minAbsValues[J][K]>curAbsValue)
						minAbsValues[J][K] = curAbsValue;
				}
			}			
		}
	}	
	
	realPrecision = pw_relBoundRatio*minAbsValues[J][K];
	floatToBytes(realPrecBytes, realPrecision);
	realPrecBytes[2] = realPrecBytes[3] = 0;
	approxPrecision = bytesToFloat(realPrecBytes);
	//put the realPrecision in float* pwrErBound
	pwrErrBound[p++] = approxPrecision;
	//put the two bytes in pwrErrBoundBytes
	pwrErrBoundBytes[q++] = realPrecBytes[0];
	pwrErrBoundBytes[q++] = realPrecBytes[1];
	
	free2DArray_float(minAbsValues, R2);
}

unsigned int optimize_intervals_float_3D_pwr(float *oriData, size_t r1, size_t r2, size_t r3, size_t R2, size_t R3, size_t edgeSize, float* pwrErrBound)
{	
	size_t i,j,k, ir,jr,index, I = 0,J=0,K=0;
	float realPrecision = pwrErrBound[0];		
	unsigned long radiusIndex;
	size_t r23=r2*r3;
	size_t R23 = R2*R3;
	float pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	size_t totalSampleSize = (r1-1)*(r2-1)*(r3-1)/sampleDistance;
	for(i=1;i<r1;i++)
	{
		ir = i*r23;
		if(i%edgeSize==0)
		{	
			I++;
			J = 0;
		}
		for(j=1;j<r2;j++)
		{
			jr = j*r3;
			if(j%edgeSize==0)
			{	
				J++;
				K = 0;
			}			
			for(k=1;k<r3;k++)
			{
				index = ir+jr+k;
				if(k%edgeSize==0)
					K++;		
				if((i+j+k)%sampleDistance==0)
				{
					realPrecision = pwrErrBound[I*R23+J*R2+K];					
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
	size_t targetCount = totalSampleSize*predThreshold;
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

void SZ_compress_args_float_NoCkRngeNoGzip_1D_pwr(unsigned char** newByteData, float *oriData, double globalPrecision, 
size_t dataLength, size_t *outSize, float min, float max)
{
	SZ_Reset();	
	size_t pwrLength = dataLength%segment_size==0?dataLength/segment_size:dataLength/segment_size+1;
	float* pwrErrBound = (float*)malloc(sizeof(float)*pwrLength);
	size_t pwrErrBoundBytes_size = sizeof(unsigned char)*pwrLength*2;
	unsigned char* pwrErrBoundBytes = (unsigned char*)malloc(pwrErrBoundBytes_size);
	
	compute_segment_precisions_float_1D(oriData, dataLength, pwrErrBound, pwrErrBoundBytes, globalPrecision);
	
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_1D_pwr(oriData, dataLength, pwrErrBound);	
		updateQuantizationInfo(quantization_intervals);
	}
	else
		quantization_intervals = intvCapacity;
	//clearHuffmanMem();
	size_t i = 0, j = 0;
	int reqLength;
	float realPrecision = pwrErrBound[j++];	
	float medianValue = 0;
	float radius = fabs(max)<fabs(min)?fabs(min):fabs(max);
	short radExpo = getExponent_float(radius);
	
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
	
	type[0] = 0;
	
	unsigned char preDataBytes[4] = {0};
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
	double predAbsErr;
	float min_pred, minErr, minIndex;
	int a = 0;
	checkRadius = (intvCapacity-1)*realPrecision;
	double interval = 2*realPrecision;
	int updateReqLength = 0; //a marker: 1 means already updated
	
	for(i=2;i<dataLength;i++)
	{
//		if(i==2499435)
//			printf("i=%d\n", i);
		curData = spaceFillingValue[i];
		if(i%segment_size==0)
		{
			realPrecision = pwrErrBound[j++];
			checkRadius = (intvCapacity-1)*realPrecision;
			interval = 2*realPrecision;
			updateReqLength = 0;
		}
		//pred = 2*last3CmprsData[0] - last3CmprsData[1];
		pred = last3CmprsData[0];
		predAbsErr = fabs(curData - pred);	
		if(predAbsErr<checkRadius)
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
		if(updateReqLength==0)
		{
			computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);
			reqBytesLength = reqLength/8;
			resiBitsLength = reqLength%8;
			updateReqLength = 1;		
		}
		
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
			realPrecision, medianValue, (char)reqLength, quantization_intervals, pwrErrBoundBytes, pwrErrBoundBytes_size, radExpo);

//sdi:Debug
/*	int sum =0;
	for(i=0;i<dataLength;i++)
		if(type[i]==0) sum++;
	printf("opt_quantizations=%d, exactDataNum=%d, sum=%d\n",quantization_intervals, exactDataNum, sum);
*/
//	writeUShortData(type, dataLength, "compressStateBytes.sb");
//	unsigned short type_[dataLength];
//	SZ_Reset();
//	decode_withTree(tdps->typeArray, tdps->typeArray_size, type_);	
//	printf("tdps->typeArray_size=%d\n", tdps->typeArray_size);
	
	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	
	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);
	
	int floatSize=sizeof(float);
	if(*outSize>dataLength*floatSize)
	{
		size_t k = 0, i;
		tdps->isLossless = 1;
		size_t totalByteLength = 3 + 4 + 1 + floatSize*dataLength;
		*newByteData = (unsigned char*)malloc(totalByteLength);
		
		unsigned char dsLengthBytes[4];
		intToBytes_bigEndian(dsLengthBytes, dataLength);//4
		for (i = 0; i < 3; i++)//3
			(*newByteData)[k++] = versionNumber[i];
		
		if(SZ_SIZE_TYPE==4)
		{
			(*newByteData)[k++] = 16;	//=00010000	
		}
		else 
		{
			(*newByteData)[k++] = 80;
		}
		for (i = 0; i < SZ_SIZE_TYPE; i++)//4 or 8
			(*newByteData)[k++] = dsLengthBytes[i];

		
		if(sysEndianType==BIG_ENDIAN_SYSTEM)
			memcpy((*newByteData)+4+SZ_SIZE_TYPE, oriData, dataLength*floatSize);
		else
		{
			unsigned char* p = (*newByteData)+4+SZ_SIZE_TYPE;
			for(i=0;i<dataLength;i++,p+=floatSize)
				floatToBytes(p, oriData[i]);
		}
		*outSize = totalByteLength;
	}

	free(pwrErrBound);
	
	free(vce);
	free(lce);
	free_TightDataPointStorageF(tdps);
	free(exactMidByteArray);
}

void SZ_compress_args_float_NoCkRngeNoGzip_2D_pwr(unsigned char** newByteData, float *oriData, double globalPrecision, size_t r1, size_t r2, 
size_t *outSize, float min, float max)
{
	SZ_Reset();	
	size_t dataLength=r1*r2;
	int blockEdgeSize = computeBlockEdgeSize_2D(segment_size);
	size_t R1 = 1+(r1-1)/blockEdgeSize;
	size_t R2 = 1+(r2-1)/blockEdgeSize;
	float* pwrErrBound = (float*)malloc(sizeof(float)*R1*R2);
	size_t pwrErrBoundBytes_size = sizeof(unsigned char)*R1*R2*2;
	unsigned char* pwrErrBoundBytes = (unsigned char*)malloc(pwrErrBoundBytes_size);
	
	compute_segment_precisions_float_2D(oriData, pwrErrBound, r1, r2, R2, blockEdgeSize, pwrErrBoundBytes, min, max, globalPrecision);
		
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{	
		quantization_intervals = optimize_intervals_float_2D_pwr(oriData, r1, r2, R2, blockEdgeSize, pwrErrBound);
		updateQuantizationInfo(quantization_intervals);
	}	
	else
		quantization_intervals = intvCapacity;
	//clearHuffmanMem();
	//printf("quantization_intervals=%d\n",quantization_intervals);
	
	size_t i=0,j=0,I=0,J=0; 
	int reqLength;
	float realPrecision = pwrErrBound[I*R2+J];	
	float pred1D, pred2D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;
	
	P0 = (float*)malloc(r2*sizeof(float));
	memset(P0, 0, r2*sizeof(float));
	P1 = (float*)malloc(r2*sizeof(float));
	memset(P1, 0, r2*sizeof(float));
		
	float medianValue = 0;
	float radius = fabs(max)<fabs(min)?fabs(min):fabs(max);	
	short radExpo = getExponent_float(radius);
	int updateReqLength = 1;
	
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
	
	type[0] = 0;
	
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
		if(j%blockEdgeSize==0)
		{
			J++;
			realPrecision = pwrErrBound[I*R2+J];
			updateReqLength = 0;
		}

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
			if(updateReqLength==0)
			{
				computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;
				updateReqLength = 1;
			}

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
	size_t index;
	for (i = 1; i < r1; i++)
	{	
		/* Process row-i data 0 */
		index = i*r2;
		J = 0;
		if(i%blockEdgeSize==0)
			I++;
		realPrecision = pwrErrBound[I*R2+J]; //J==0
		updateReqLength = 0;
		
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
			if(updateReqLength==0)
			{
				computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;
				updateReqLength = 1;
			}
			
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
			
			if(j%blockEdgeSize==0)
			{
				J++;
				realPrecision = pwrErrBound[I*R2+J];
				updateReqLength = 0;
			}
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
				if(updateReqLength==0)
				{
					computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);
					reqBytesLength = reqLength/8;
					resiBitsLength = reqLength%8;
					updateReqLength = 1;
				}

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
			realPrecision, medianValue, (char)reqLength, quantization_intervals, pwrErrBoundBytes, pwrErrBoundBytes_size, radExpo);
	
	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);
	
	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);
	
	free(pwrErrBound);

	free(vce);
	free(lce);
	free_TightDataPointStorageF(tdps);	
	free(exactMidByteArray);
}

void SZ_compress_args_float_NoCkRngeNoGzip_3D_pwr(unsigned char** newByteData, float *oriData, double globalPrecision, 
size_t r1, size_t r2, size_t r3, size_t *outSize, float min, float max)
{
	SZ_Reset();	
	size_t dataLength=r1*r2*r3;
	
	int blockEdgeSize = computeBlockEdgeSize_3D(segment_size);
	size_t R1 = 1+(r1-1)/blockEdgeSize;
	size_t R2 = 1+(r2-1)/blockEdgeSize;
	size_t R3 = 1+(r3-1)/blockEdgeSize;
	//float*** pwrErrBound = create3DArray_float(R1, R2, R3);
	float* pwrErrBound = (float*)malloc(sizeof(float)*R1*R2*R3);
	size_t pwrErrBoundBytes_size = sizeof(unsigned char)*R1*R2*R3*2;
	unsigned char* pwrErrBoundBytes = (unsigned char*)malloc(pwrErrBoundBytes_size);	
	
	compute_segment_precisions_float_3D(oriData, pwrErrBound, r1, r2, r3, R2, R3, blockEdgeSize, pwrErrBoundBytes, min, max, globalPrecision);	

	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_float_3D_pwr(oriData, r1, r2, r3, R2, R3, blockEdgeSize, pwrErrBound);
		updateQuantizationInfo(quantization_intervals);
	}	
	else
		quantization_intervals = intvCapacity;
	//clearHuffmanMem();
	size_t i=0,j=0,k=0, I = 0, J = 0, K = 0;
	int reqLength;
	float realPrecision = pwrErrBound[0];		
	float pred1D, pred2D, pred3D;
	float diff = 0.0;
	double itvNum = 0;
	float *P0, *P1;

	size_t r23 = r2*r3;
	size_t R23 = R2*R3;
	P0 = (float*)malloc(r23*sizeof(float));
	P1 = (float*)malloc(r23*sizeof(float));
	float radius = fabs(max)<fabs(min)?fabs(min):fabs(max);
	float medianValue = 0;
	short radExpo = getExponent_float(radius);
	int updateReqLength = 0;
	
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

	type[0] = 0;

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
		if(updateReqLength==0)
		{
			computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);
			reqBytesLength = reqLength/8;
			resiBitsLength = reqLength%8;
			updateReqLength = 1;
		}		
		
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
		if(j%blockEdgeSize==0)
		{
			J++;
			realPrecision = pwrErrBound[J];
			updateReqLength = 0;
		}		
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
			if(updateReqLength==0)
			{
				computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;
				updateReqLength = 1;
			}			

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
	size_t index;
	K = 0;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;	

		J = 0;
		if(i%blockEdgeSize==0)
			I++;
		realPrecision = pwrErrBound[I*R3+J]; //J==0
		updateReqLength = 0;

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
			if(updateReqLength==0)
			{
				computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;
				updateReqLength = 1;
			}		
						
			type[index] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleFloatValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index] = vce->data;
		}

		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++) //note that this j refers to fastest dimension (lowest order)
		{
			index = i*r3+j;		
			if(j%blockEdgeSize==0)
			{
				J++;
				realPrecision = pwrErrBound[I*R3+J];
				updateReqLength = 0;
			}			
		
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
				if(updateReqLength==0)
				{
					computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);
					reqBytesLength = reqLength/8;
					resiBitsLength = reqLength%8;
					updateReqLength = 1;
				}						
				
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
		I = 0;
		J = 0;
		if(k%blockEdgeSize==0)
			K++;
		realPrecision = pwrErrBound[K*R23]; //J==0
		updateReqLength = 0;
		
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
			if(updateReqLength==0)
			{
				computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;
				updateReqLength = 1;
			}					
			
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
			index = k*r23+j;	

			if(j%blockEdgeSize==0)
			{
				J++;
				realPrecision = pwrErrBound[K*R23+J];
				updateReqLength = 0;			
			}					
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
				if(updateReqLength==0)
				{
					computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);
					reqBytesLength = reqLength/8;
					resiBitsLength = reqLength%8;
					updateReqLength = 1;
				}						
				
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
		size_t index2D;
		for (i = 1; i < r2; i++)
		{
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;
			J = 0;
			if(i%blockEdgeSize==0)
				I++;
			realPrecision = pwrErrBound[K*R23+I*R3+J]; //J==0
			updateReqLength = 0;			
			
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
				if(updateReqLength==0)
				{
					computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);
					reqBytesLength = reqLength/8;
					resiBitsLength = reqLength%8;
					updateReqLength = 1;
				}						
				
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
				index = k*r23 + i*r3 + j;
				if(j%blockEdgeSize==0)
				{
					J++;
					realPrecision = pwrErrBound[K*R23+I*R3+J];
					updateReqLength = 0;			
				}							
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
					if(updateReqLength==0)
					{
						computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);
						reqBytesLength = reqLength/8;
						resiBitsLength = reqLength%8;
						updateReqLength = 1;
					}							
					
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
			realPrecision, medianValue, (char)reqLength, quantization_intervals, pwrErrBoundBytes, pwrErrBoundBytes_size, radExpo);

//sdi:Debug
/*	int sum =0;
	for(i=0;i<dataLength;i++)
		if(type[i]==0) sum++;
	printf("opt_quantizations=%d, exactDataNum=%d, sum=%d\n",quantization_intervals, exactDataNum, sum);
*/

	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);

	convertTDPStoFlatBytes_float(tdps, newByteData, outSize);

	free(pwrErrBound);

	free(vce);
	free(lce);
	free_TightDataPointStorageF(tdps);
	free(exactMidByteArray);
}

void createRangeGroups_float(float** posGroups, float** negGroups, int** posFlags, int** negFlags)
{
	size_t size = GROUP_COUNT*sizeof(float);
	size_t size2 = GROUP_COUNT*sizeof(int);
	*posGroups = (float*)malloc(size);
	*negGroups = (float*)malloc(size);
	*posFlags = (int*)malloc(size2);
	*negFlags = (int*)malloc(size2);
	memset(*posGroups, 0, size);
	memset(*negGroups, 0, size);
	memset(*posFlags, 0, size2);
	memset(*negFlags, 0, size2);
}

void compressGroupIDArray_float(char* groupID, TightDataPointStorageF* tdps)
{
	SZ_Reset(); //reset the huffman tree
	size_t dataLength = tdps->dataSeriesLength;
	int* standGroupID = (int*)malloc(dataLength*sizeof(int));

	size_t i;
	standGroupID[0] = groupID[0]+GROUP_COUNT; //plus an offset such that it would not be a negative number.
	char lastGroupIDValue = groupID[0], curGroupIDValue;
	int offset = 2*(GROUP_COUNT + 2);
	for(i=1; i<dataLength;i++)
	{
		curGroupIDValue = groupID[i];
		standGroupID[i] = (curGroupIDValue - lastGroupIDValue) + offset; 
		lastGroupIDValue = curGroupIDValue;
	}
	
	unsigned char* out = (unsigned char*)malloc(sizeof(unsigned char)*dataLength);
	size_t outSize;
	
	encode_withTree(standGroupID, dataLength, &out, &outSize);
	//encode2_withTree(standGroupID, dataLength, &out, &outSize);
	
	//outSize = zlib_compress5(standGroupID, dataLength, &out, 9);
	
	tdps->pwrErrBoundBytes = out; //groupIDArray
	tdps->pwrErrBoundBytes_size = outSize;
	
	free(standGroupID);
}

TightDataPointStorageF* SZ_compress_float_1D_MDQ_pwrGroup(float* oriData, size_t dataLength, int errBoundMode, 
double absErrBound, double relBoundRatio, double pwrErrRatio, float valueRangeSize, float medianValue_f)
{
	size_t i;
	float *posGroups, *negGroups, *groups;
	float pos_01_group = 0, neg_01_group = 0; //[0,1] and [-1,0]
	int *posFlags, *negFlags, *flags;
	int pos_01_flag = 0, neg_01_flag = 0;
	createRangeGroups_float(&posGroups, &negGroups, &posFlags, &negFlags);
	size_t nbBins = (size_t)(1/pwrErrRatio);
	if(nbBins%2==1)
		nbBins++;
	intvRadius = nbBins;

	int reqLength, status;
	float medianValue = medianValue_f;
	float realPrecision = (float)getRealPrecision_float(valueRangeSize, errBoundMode, absErrBound, relBoundRatio, &status);
	if(realPrecision<0)
		realPrecision = pwrErrRatio;
	float realGroupPrecision; //precision (error) based on group ID
	short reqExpo = getPrecisionReqLength_float(realPrecision);
	short radExpo = getExponent_float(valueRangeSize/2);
	short lastGroupNum, groupNum, decGroupNum, grpNum = 0;
	
	double* groupErrorBounds = generateGroupErrBounds(errBoundMode, realPrecision, pwrErrRatio);
	intvRadius = generateGroupMaxIntervalCount(groupErrorBounds);
	
	computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));
	char *groupID = (char*) malloc(dataLength*sizeof(char));
	char *gp = groupID;
		
	float* spaceFillingValue = oriData; 
	
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
			
	int state;
	float lcf, qcf;	
	float curData, decValue;
	float pred;
	float predAbsErr;
	double errBound, interval = 0;
	
	//add the first data	
	type[0] = 0;
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
	compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	
	curData = spaceFillingValue[0];
	groupNum = computeGroupNum_float(vce->data);

	if(curData > 0 && groupNum >= 0)
	{
		groups = posGroups;
		flags = posFlags;
		grpNum = groupNum;
	}
	else if(curData < 0 && groupNum >= 0)
	{
		groups = negGroups;
		flags = negFlags;
		grpNum = groupNum;
	}
	else if(curData >= 0 && groupNum == -1)
	{
		groups = &pos_01_group;
		flags = &pos_01_flag;
		grpNum = 0;
	}
	else //curData < 0 && groupNum == -1
	{
		groups = &neg_01_group;
		flags = &neg_01_flag;
		grpNum = 0;
	}

	listAdd_float_group(groups, flags, groupNum, spaceFillingValue[0], vce->data, gp);
	gp++;
	
	for(i=1;i<dataLength;i++)
	{
		curData = oriData[i];
		//printf("i=%d, posGroups[3]=%f, negGroups[3]=%f\n", i, posGroups[3], negGroups[3]);
		
		groupNum = computeGroupNum_float(curData);
		
		if(curData > 0 && groupNum >= 0)
		{
			groups = posGroups;
			flags = posFlags;
			grpNum = groupNum;
		}
		else if(curData < 0 && groupNum >= 0)
		{
			groups = negGroups;
			flags = negFlags;
			grpNum = groupNum;
		}
		else if(curData >= 0 && groupNum == -1)
		{
			groups = &pos_01_group;
			flags = &pos_01_flag;
			grpNum = 0;
		}
		else //curData < 0 && groupNum == -1
		{
			groups = &neg_01_group;
			flags = &neg_01_flag;
			grpNum = 0;
		}

		if(groupNum>=GROUP_COUNT)
		{
			type[i] = 0;
			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);

			compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			listAdd_float_group(groups, flags, lastGroupNum, curData, vce->data, gp);	//set the group number to be last one in order to get the groupID array as smooth as possible.		
		}
		else if(flags[grpNum]==0) //the dec value may not be in the same group
		{	
			type[i] = 0;
			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);

			compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,4);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			//decGroupNum = computeGroupNum_float(vce->data);
			
			//if(decGroupNum < groupNum)
			//	decValue = curData>0?pow(2, groupNum):-pow(2, groupNum);
			//else if(decGroupNum > groupNum)
			//	decValue = curData>0?pow(2, groupNum+1):-pow(2, groupNum+1);
			//else
			//	decValue = vce->data;
			
			decValue = vce->data;	
			listAdd_float_group(groups, flags, groupNum, curData, decValue, gp);
			lastGroupNum = curData>0?groupNum + 2: -(groupNum+2);
		}
		else //if flags[groupNum]==1, the dec value must be in the same group
		{
			pred = groups[grpNum];
			predAbsErr = fabs(curData - pred);
			realGroupPrecision = groupErrorBounds[grpNum]; //compute real error bound
			interval = realGroupPrecision*2;
			state = (predAbsErr/realGroupPrecision+1)/2;
			if(curData>=pred)
			{
				type[i] = intvRadius+state;
				decValue = pred + state*interval;
			}
			else //curData<pred
			{
				type[i] = intvRadius-state;
				decValue = pred - state*interval;
			}
			//decGroupNum = computeGroupNum_float(pred);
			
			if((decValue>0&&curData<0)||(decValue<0&&curData>=0))
				decValue = 0;
			//else
			//{
			//	if(decGroupNum < groupNum)
			//		decValue = curData>0?pow(2, groupNum):-pow(2, groupNum);
			//	else if(decGroupNum > groupNum)
			//		decValue = curData>0?pow(2, groupNum+1):-pow(2, groupNum+1);
			//	else
			//		decValue = pred;				
			//}
			
			if(fabs(curData-decValue)>realGroupPrecision)
			{	
				type[i] = 0;
				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);

				compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,4);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);

				decValue = vce->data;	
			}
			
			listAdd_float_group(groups, flags, groupNum, curData, decValue, gp);			
			lastGroupNum = curData>=0?groupNum + 2: -(groupNum+2);			
		}
		gp++;	

	}
	
	int exactDataNum = exactLeadNumArray->size;
	
	TightDataPointStorageF* tdps;
			
	//combineTypeAndGroupIDArray(nbBins, dataLength, &type, groupID);

	new_TightDataPointStorageF(&tdps, dataLength, exactDataNum, 
			type, exactMidByteArray->array, exactMidByteArray->size,  
			exactLeadNumArray->array,  
			resiBitArray->array, resiBitArray->size, 
			resiBitLengthArray->array, resiBitLengthArray->size, 
			realPrecision, medianValue, (char)reqLength, nbBins, NULL, 0, radExpo);	
	
	compressGroupIDArray_float(groupID, tdps);
	
	free(posGroups);
	free(negGroups);
	free(posFlags);
	free(negFlags);
	free(groupID);
	free(groupErrorBounds);
	
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);	
	free(vce);
	free(lce);	
	free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);	
	
	return tdps;
}

void SZ_compress_args_float_NoCkRngeNoGzip_1D_pwrgroup(unsigned char** newByteData, float *oriData,
size_t dataLength, double absErrBound, double relBoundRatio, double pwrErrRatio, float valueRangeSize, float medianValue_f, size_t *outSize)
{
        SZ_Reset();
        TightDataPointStorageF* tdps = SZ_compress_float_1D_MDQ_pwrGroup(oriData, dataLength, errorBoundMode, 
        absErrBound, relBoundRatio, pwrErrRatio, 
        valueRangeSize, medianValue_f);

        convertTDPStoFlatBytes_float(tdps, newByteData, outSize);

        if(*outSize>dataLength*sizeof(float))
                SZ_compress_args_float_StoreOriData(oriData, dataLength+2, tdps, newByteData, outSize);

        free_TightDataPointStorageF(tdps);
}

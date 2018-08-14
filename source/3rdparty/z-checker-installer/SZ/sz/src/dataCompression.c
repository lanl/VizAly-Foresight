/**
 *  @file double_compression.c
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Compression Technique for double array
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "sz.h"
#include "DynamicByteArray.h"
#include "DynamicIntArray.h"
#include "TightDataPointStorageD.h"
#include "CompressElement.h"

float computeRangeSize_float(float* oriData, int size, float* valueRangeSize, float* medianValue)
{
	int i = 0;
	float min = oriData[0];
	float max = min;
	for(i=1;i<size;i++)
	{
		float data = oriData[i];
		if(min>data)
			min = data;
		else if(max<data)
			max = data;
	}

	*valueRangeSize = max - min;
	*medianValue = min + *valueRangeSize/2;
	return min;
}

double computeRangeSize_double(double* oriData, int size, double* valueRangeSize, double* medianValue)
{
	int i = 0;
	double min = oriData[0];
	double max = min;
	for(i=1;i<size;i++)
	{
		double data = oriData[i];
		if(min>data)
			min = data;
		else if(max<data)
			max = data;
	}
	
	*valueRangeSize = max - min;
	*medianValue = min + *valueRangeSize/2;
	return min;
}

float computeRangeSize_float_subblock(float* oriData, float* valueRangeSize, float* medianValue,
int r5, int r4, int r3, int r2, int r1,
int s5, int s4, int s3, int s2, int s1,
int e5, int e4, int e3, int e2, int e1)
{
	int i1, i2, i3, i4, i5;
	int index_start = s5*(r4*r3*r2*r1) + s4*(r3*r2*r1) + s3*(r2*r1) + s2*r1 + s1;
	float min = oriData[index_start];
	float max = min;

	for (i5 = s5; i5 <= e5; i5++)
	for (i4 = s4; i4 <= e4; i4++)
	for (i3 = s3; i3 <= e3; i3++)
	for (i2 = s2; i2 <= e2; i2++)
	for (i1 = s1; i1 <= e1; i1++)
	{
		int index = i5*(r4*r3*r2*r1) + i4*(r3*r2*r1) + i3*(r2*r1) + i2*r1 + i1;
		float data = oriData[index];
		if (min>data)
			min = data;
		else if(max<data)
			max = data;
	}

	*valueRangeSize = max - min;
	*medianValue = min + *valueRangeSize/2;
	return min;
}


float computeRangeSize_double_subblock(double* oriData, double* valueRangeSize, double* medianValue,
int r5, int r4, int r3, int r2, int r1,
int s5, int s4, int s3, int s2, int s1,
int e5, int e4, int e3, int e2, int e1)
{
	int i1, i2, i3, i4, i5;
	int index_start = s5*(r4*r3*r2*r1) + s4*(r3*r2*r1) + s3*(r2*r1) + s2*r1 + s1;
	double min = oriData[index_start];
	double max = min;

	for (i5 = s5; i5 <= e5; i5++)
	for (i4 = s4; i4 <= e4; i4++)
	for (i3 = s3; i3 <= e3; i3++)
	for (i2 = s2; i2 <= e2; i2++)
	for (i1 = s1; i1 <= e1; i1++)
	{
		int index = i5*(r4*r3*r2*r1) + i4*(r3*r2*r1) + i3*(r2*r1) + i2*r1 + i1;
		double data = oriData[index];
		if (min>data)
			min = data;
		else if(max<data)
			max = data;
	}

	*valueRangeSize = max - min;
	*medianValue = min + *valueRangeSize/2;
	return min;
}


double min_d(double a, double b)
{
	if(a<b)
		return a;
	else
		return b;
}

double max_d(double a, double b)
{
	if(a>b)
		return a;
	else
		return b;
}

float min_f(float a, float b)
{
	if(a<b)
		return a;
	else
		return b;
}

float max_f(float a, float b)
{
	if(a>b)
		return a;
	else
		return b;
}

double getRealPrecision_double(double valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status)
{
	int state = SZ_SCES;
	double precision = 0;
	if(errBoundMode==ABS)
		precision = absErrBound; 
	else if(errBoundMode==REL)
		precision = relBoundRatio*valueRangeSize;
	else if(errBoundMode==ABS_AND_REL)
		precision = min_d(absErrBound, relBoundRatio*valueRangeSize);
	else if(errBoundMode==ABS_OR_REL)
		precision = max_d(absErrBound, relBoundRatio*valueRangeSize);
	else if(errBoundMode==PW_REL)
		precision = -1;
	else
	{
		printf("Error: error-bound-mode is incorrect!\n");
		state = SZ_BERR;
	}
	*status = state;
	return precision;
}

double getRealPrecision_float(float valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status)
{
	int state = SZ_SCES;
	double precision = 0;
	if(errBoundMode==ABS)
		precision = absErrBound; 
	else if(errBoundMode==REL)
		precision = relBoundRatio*valueRangeSize;
	else if(errBoundMode==ABS_AND_REL)
		precision = min_f(absErrBound, relBoundRatio*valueRangeSize);
	else if(errBoundMode==ABS_OR_REL)
		precision = max_f(absErrBound, relBoundRatio*valueRangeSize);
	else if(errBoundMode==PW_REL)
		precision = -1;
	else
	{
		printf("Error: error-bound-mode is incorrect!\n");
		state = SZ_BERR;
	}
	*status = state;
	return precision;
}

void symTransform_8bytes(unsigned char data[8])
{
	unsigned char tmp = data[0];
	data[0] = data[7];
	data[7] = tmp;

	tmp = data[1];
	data[1] = data[6];
	data[6] = tmp;
	
	tmp = data[2];
	data[2] = data[5];
	data[5] = tmp;
	
	tmp = data[3];
	data[3] = data[4];
	data[4] = tmp;
}

void flush_to_bigEndian_8bytes(unsigned char data[8], int dataEndianType)
{
	if(dataEndianType==LITTLE_ENDIAN_DATA)
		symTransform_8bytes(data);
}

inline void symTransform_2bytes(unsigned char data[2])
{
	unsigned char tmp = data[0];
	data[0] = data[1];
	data[1] = tmp;
}

void symTransform_4bytes(unsigned char data[4])
{
	unsigned char tmp = data[0];
	data[0] = data[3];
	data[3] = tmp;

	tmp = data[1];
	data[1] = data[2];
	data[2] = tmp;
}

void flush_to_bigEndian_4bytes(unsigned char data[4])
{
	if(dataEndianType==LITTLE_ENDIAN_DATA)
		symTransform_4bytes(data);
}

void bigEndian_to_OSEndian_double(unsigned char data[8])
{
	if(dataEndianType==LITTLE_ENDIAN_DATA)
		symTransform_8bytes(data);
}

void bigEndian_to_OSEndian_float(unsigned char data[4])
{
	if(dataEndianType==LITTLE_ENDIAN_DATA)
		symTransform_4bytes(data);
}

void compressSingleFloatValue(FloatValueCompressElement *vce, float tgtValue, float precision, float medianValue, 
		int reqLength, int reqBytesLength, int resiBitsLength)
{		
	float normValue = tgtValue - medianValue;

	lfloat lfBuf;
	lfBuf.value = normValue;
			
	int ignBytesLength = 32 - reqLength;
	if(ignBytesLength<0)
		ignBytesLength = 0;
	
	int tmp_int = lfBuf.ivalue;
	intToBytes_bigEndian(vce->curBytes, tmp_int);
		
	lfBuf.ivalue = (lfBuf.ivalue >> ignBytesLength) << ignBytesLength;
	
	//float tmpValue = lfBuf.value;
	
	vce->data = lfBuf.value+medianValue;
	vce->curValue = tmp_int;
	vce->reqBytesLength = reqBytesLength;
	vce->resiBitsLength = resiBitsLength;
}

void compressSingleDoubleValue(DoubleValueCompressElement *vce, double tgtValue, double precision, double medianValue, 
		int reqLength, int reqBytesLength, int resiBitsLength)
{		
	double normValue = tgtValue - medianValue;

	ldouble lfBuf;
	lfBuf.value = normValue;
			
	int ignBytesLength = 64 - reqLength;
	if(ignBytesLength<0)
		ignBytesLength = 0;

	long tmp_long = lfBuf.lvalue;
	longToBytes_bigEndian(vce->curBytes, tmp_long);
				
	lfBuf.lvalue = (lfBuf.lvalue >> ignBytesLength)<<ignBytesLength;
	
	//double tmpValue = lfBuf.value;
	
	vce->data = lfBuf.value+medianValue;
	vce->curValue = tmp_long;
	vce->reqBytesLength = reqBytesLength;
	vce->resiBitsLength = resiBitsLength;
}

int compIdenticalLeadingBytesCount_double(unsigned char* preBytes, unsigned char* curBytes)
{
	int i, n = 0;
	for(i=0;i<8;i++)
		if(preBytes[i]==curBytes[i])
			n++;
		else
			break;
	if(n>3) n = 3;
	return n;
}

int compIdenticalLeadingBytesCount_float(unsigned char* preBytes, unsigned char* curBytes)
{
	int i, n = 0;
	for(i=0;i<4;i++)
		if(preBytes[i]==curBytes[i])
			n++;
		else
			break;
	if(n>3) n = 3;
	return n;
}

//TODO double-check the correctness...
void addExactData(DynamicByteArray *exactMidByteArray, DynamicIntArray *exactLeadNumArray, 
		DynamicIntArray *resiBitArray, LossyCompressionElement *lce)
{
	int i;
	int leadByteLength = lce->leadingZeroBytes;
	addDIA_Data(exactLeadNumArray, leadByteLength);
	unsigned char* intMidBytes = lce->integerMidBytes;
	int integerMidBytesLength = lce->integerMidBytes_Length;
	int resMidBitsLength = lce->resMidBitsLength;
	if(intMidBytes!=NULL||resMidBitsLength!=0)
	{
		if(intMidBytes!=NULL)
			for(i = 0;i<integerMidBytesLength;i++)
				addDBA_Data(exactMidByteArray, intMidBytes[i]);
		if(resMidBitsLength!=0)
			addDIA_Data(resiBitArray, lce->residualMidBits);
	}
}

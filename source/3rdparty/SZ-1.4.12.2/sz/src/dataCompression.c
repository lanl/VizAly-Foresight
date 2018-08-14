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
#include "dataCompression.h"

int computeByteSizePerIntValue(long valueRangeSize)
{
	if(valueRangeSize<=256)
		return 1;
	else if(valueRangeSize<=65536)
		return 2;
	else if(valueRangeSize<=4294967296) //2^32
		return 4;
	else
		return 8;
}

long computeRangeSize_int(void* oriData, int dataType, size_t size, int64_t* valueRangeSize)
{
	size_t i = 0;
	long max, min;

	if(dataType==SZ_UINT8)
	{
		unsigned char* data = (unsigned char*)oriData;
		unsigned char data_; 
		min = data[0], max = min;
		computeMinMax(data);
	}
	else if(dataType == SZ_INT8)
	{
		char* data = (char*)oriData;
		char data_;
		min = data[0], max = min;
		computeMinMax(data);
	}
	else if(dataType == SZ_UINT16)
	{
		unsigned short* data = (unsigned short*)oriData;
		unsigned short data_; 
		min = data[0], max = min;
		computeMinMax(data);
	}
	else if(dataType == SZ_INT16)
	{ 
		short* data = (short*)oriData;
		short data_; 
		min = data[0], max = min;
		computeMinMax(data);
	}
	else if(dataType == SZ_UINT32)
	{
		unsigned int* data = (unsigned int*)oriData;
		int data_; 
		min = data[0], max = min;
		computeMinMax(data);
	}
	else if(dataType == SZ_INT32)
	{
		int* data = (int*)oriData;
		unsigned int data_; 
		min = data[0], max = min;
		computeMinMax(data);
	}
	else if(dataType == SZ_UINT64)
	{
		unsigned long* data = (unsigned long*)oriData;
		unsigned long data_; 
		min = data[0], max = min;
		computeMinMax(data);
	}
	else if(dataType == SZ_INT64)
	{
		long* data = (long *)oriData;
		long data_; 
		min = data[0], max = min;
		computeMinMax(data);
	}

	*valueRangeSize = max - min;
	return min;	
}

float computeRangeSize_float(float* oriData, size_t size, float* valueRangeSize, float* medianValue)
{
	size_t i = 0;
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

double computeRangeSize_double(double* oriData, size_t size, double* valueRangeSize, double* medianValue)
{
	size_t i = 0;
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
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
size_t s5, size_t s4, size_t s3, size_t s2, size_t s1,
size_t e5, size_t e4, size_t e3, size_t e2, size_t e1)
{
	size_t i1, i2, i3, i4, i5;
	size_t index_start = s5*(r4*r3*r2*r1) + s4*(r3*r2*r1) + s3*(r2*r1) + s2*r1 + s1;
	float min = oriData[index_start];
	float max = min;

	for (i5 = s5; i5 <= e5; i5++)
	for (i4 = s4; i4 <= e4; i4++)
	for (i3 = s3; i3 <= e3; i3++)
	for (i2 = s2; i2 <= e2; i2++)
	for (i1 = s1; i1 <= e1; i1++)
	{
		size_t index = i5*(r4*r3*r2*r1) + i4*(r3*r2*r1) + i3*(r2*r1) + i2*r1 + i1;
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
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
size_t s5, size_t s4, size_t s3, size_t s2, size_t s1,
size_t e5, size_t e4, size_t e3, size_t e2, size_t e1)
{
	size_t i1, i2, i3, i4, i5;
	size_t index_start = s5*(r4*r3*r2*r1) + s4*(r3*r2*r1) + s3*(r2*r1) + s2*r1 + s1;
	double min = oriData[index_start];
	double max = min;

	for (i5 = s5; i5 <= e5; i5++)
	for (i4 = s4; i4 <= e4; i4++)
	for (i3 = s3; i3 <= e3; i3++)
	for (i2 = s2; i2 <= e2; i2++)
	for (i1 = s1; i1 <= e1; i1++)
	{
		size_t index = i5*(r4*r3*r2*r1) + i4*(r3*r2*r1) + i3*(r2*r1) + i2*r1 + i1;
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
	if(errBoundMode==ABS||errBoundMode==ABS_OR_PW_REL||errBoundMode==ABS_AND_PW_REL)
		precision = absErrBound; 
	else if(errBoundMode==REL||errBoundMode==REL_OR_PW_REL||errBoundMode==REL_AND_PW_REL)
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
	if(errBoundMode==ABS||errBoundMode==ABS_OR_PW_REL||errBoundMode==ABS_AND_PW_REL)
		precision = absErrBound; 
	else if(errBoundMode==REL||errBoundMode==REL_OR_PW_REL||errBoundMode==REL_AND_PW_REL)
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

double getRealPrecision_int(long valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status)
{
	int state = SZ_SCES;
	double precision = 0;
	if(errBoundMode==ABS||errBoundMode==ABS_OR_PW_REL||errBoundMode==ABS_AND_PW_REL)
		precision = absErrBound; 
	else if(errBoundMode==REL||errBoundMode==REL_OR_PW_REL||errBoundMode==REL_AND_PW_REL)
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

inline void symTransform_2bytes(unsigned char data[2])
{
	unsigned char tmp = data[0];
	data[0] = data[1];
	data[1] = tmp;
}

inline void symTransform_4bytes(unsigned char data[4])
{
	unsigned char tmp = data[0];
	data[0] = data[3];
	data[3] = tmp;

	tmp = data[1];
	data[1] = data[2];
	data[2] = tmp;
}

inline void compressInt8Value(int8_t tgtValue, int8_t minValue, int byteSize, unsigned char* bytes)
{
	uint8_t data = tgtValue - minValue;
	memcpy(bytes, &data, byteSize); //byteSize==1
}

inline void compressInt16Value(int16_t tgtValue, int16_t minValue, int byteSize, unsigned char* bytes)
{
	uint16_t data = tgtValue - minValue;
	unsigned char tmpBytes[2];
	int16ToBytes_bigEndian(tmpBytes, data);
	memcpy(bytes, tmpBytes + 2 - byteSize, byteSize);
}

inline void compressInt32Value(int32_t tgtValue, int32_t minValue, int byteSize, unsigned char* bytes)
{
	uint32_t data = tgtValue - minValue;
	unsigned char tmpBytes[4];
	int32ToBytes_bigEndian(tmpBytes, data);
	memcpy(bytes, tmpBytes + 4 - byteSize, byteSize);
}

inline void compressInt64Value(int64_t tgtValue, int64_t minValue, int byteSize, unsigned char* bytes)
{
	uint64_t data = tgtValue - minValue;
	unsigned char tmpBytes[8];
	int64ToBytes_bigEndian(tmpBytes, data);
	memcpy(bytes, tmpBytes + 8 - byteSize, byteSize);
}

inline void compressUInt8Value(uint8_t tgtValue, uint8_t minValue, int byteSize, unsigned char* bytes)
{
	uint8_t data = tgtValue - minValue;
	memcpy(bytes, &data, byteSize); //byteSize==1
}

inline void compressUInt16Value(uint16_t tgtValue, uint16_t minValue, int byteSize, unsigned char* bytes)
{
	uint16_t data = tgtValue - minValue;
	unsigned char tmpBytes[2];
	int16ToBytes_bigEndian(tmpBytes, data);
	memcpy(bytes, tmpBytes + 2 - byteSize, byteSize);
}

inline void compressUInt32Value(uint32_t tgtValue, uint32_t minValue, int byteSize, unsigned char* bytes)
{
	uint32_t data = tgtValue - minValue;
	unsigned char tmpBytes[4];
	int32ToBytes_bigEndian(tmpBytes, data);
	memcpy(bytes, tmpBytes + 4 - byteSize, byteSize);
}

inline void compressUInt64Value(uint64_t tgtValue, uint64_t minValue, int byteSize, unsigned char* bytes)
{
	uint64_t data = tgtValue - minValue;
	unsigned char tmpBytes[8];
	int64ToBytes_bigEndian(tmpBytes, data);
	memcpy(bytes, tmpBytes + 8 - byteSize, byteSize);
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

/**
 * @deprecated
 * @return: the length of the coefficient array.
 * */
int getPredictionCoefficients(int layers, int dimension, int **coeff_array, int *status)
{
	size_t size = 0;
	switch(dimension)
	{
		case 1:
			switch(layers)
			{
				case 1:
					*coeff_array = (int*)malloc(sizeof(int));
					(*coeff_array)[0] = 1;
					size = 1;
					break;
				case 2:
					*coeff_array = (int*)malloc(2*sizeof(int));
					(*coeff_array)[0] = 2;
					(*coeff_array)[1] = -1;
					size = 2;
					break;
				case 3:
					*coeff_array = (int*)malloc(3*sizeof(int));
					(*coeff_array)[0] = 3;
					(*coeff_array)[1] = -3;
					(*coeff_array)[2] = 1;
					break;
			}	
			break;
		case 2:
			switch(layers)
			{
				case 1:
				
					break;
				case 2:
				
					break;
				case 3:
				
					break;
			}				
			break;
		case 3:
			switch(layers)
			{
				case 1:
				
					break;
				case 2:
				
					break;
				case 3:
				
					break;
			}			
			break;
		default:
			printf("Error: dimension must be no greater than 3 in the current version.\n");
			*status = SZ_DERR;
	}
	*status = SZ_SCES;
	return size;
}

int computeBlockEdgeSize_2D(int segmentSize)
{
	int i = 1;
	for(i=1; i<segmentSize;i++)
	{
		if(i*i>segmentSize)
			break;
	}
	return i;
	//return (int)(sqrt(segmentSize)+1);
}

int computeBlockEdgeSize_3D(int segmentSize)
{
	int i = 1;
	for(i=1; i<segmentSize;i++)
	{
		if(i*i*i>segmentSize)
			break;
	}
	return i;	
	//return (int)(pow(segmentSize, 1.0/3)+1);
}

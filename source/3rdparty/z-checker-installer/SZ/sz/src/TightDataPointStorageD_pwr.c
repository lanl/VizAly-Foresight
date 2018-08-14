/**
 *  @file TightPointDataStorageD.c
 *  @author Sheng Di
 *  @date May, 2016
 *  @brief The functions used to construct the tightPointDataStorage element for storing compressed bytes.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "TightDataPointStorageD.h"
#include "sz.h"
#include "Huffman.h"
//#include "rw.h"

void decompressDataSeries_double_1D_pwr(double** data, int dataSeriesLength, TightDataPointStorageD* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	unsigned char tmpPrecBytes[8] = {0}; //used when needing to convert bytes to float values
	unsigned char* bp = tdps->pwrErrBoundBytes;
	int i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
								// in resiMidBits, p is to track the
								// byte_index of resiMidBits, l is for
								// leadNum

	unsigned char* leadNum;
	double interval;// = (float)tdps->realPrecision*2;
	
	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);
	*data = (double*)malloc(sizeof(double)*dataSeriesLength);

	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);
	
	unsigned char preBytes[8];
	unsigned char curBytes[8];
	
	memset(preBytes, 0, 8);

	int curByteIndex = 0;
	int reqLength, reqBytesLength, resiBitsLength, resiBits, reqExpo, reqMantLength; 
	unsigned char leadingNum;	
	double medianValue, exactData, predValue, realPrecision;
	
	medianValue = tdps->medianValue;
	
	int type_, updateReqLength = 0;
	for (i = 0; i < dataSeriesLength; i++) 
	{
		if(i%tdps->segment_size==0)
		{
			tmpPrecBytes[0] = *(bp++);
			tmpPrecBytes[1] = *(bp++);
			memset(&tmpPrecBytes[2], 0, 6*sizeof(unsigned char));

			realPrecision = bytesToDouble(tmpPrecBytes);
			interval = realPrecision*2;
			updateReqLength = 0;
		}
		type_ = type[i];
		switch (type_) {
		case 0:
			// compute resiBits
			if(updateReqLength==0)
			{
				computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;	
				updateReqLength = 1;	
			}
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}
			
			exactData = bytesToDouble(curBytes);
			(*data)[i] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
			break;
		default:
			//predValue = 2 * (*data)[i-1] - (*data)[i-2];
			predValue = (*data)[i-1];
			(*data)[i] = predValue + (type_-intvRadius)*interval;
			break;
		}
		//printf("%.30G\n",(*data)[i]);
	}
	free(leadNum);
	free(type);
	return;
}

double* extractRealPrecision_2D_double(int R1, int R2, int blockSize, TightDataPointStorageD* tdps)
{
	int i,j,k=0, I;
	unsigned char* bytes = tdps->pwrErrBoundBytes;
	unsigned char tmpBytes[8] = {0};
	double* result = (double*)malloc(sizeof(double)*R1*R2);
	for(i=0;i<R1;i++)
	{
		I = i*R2;
		for(j=0;j<R2;j++)
		{
			tmpBytes[0] = bytes[k++];
			tmpBytes[1] = bytes[k++];
			result[I+j]=bytesToDouble(tmpBytes);
		}
	}
	return result;
}

void decompressDataSeries_double_2D_pwr(double** data, int r1, int r2, TightDataPointStorageD* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	//printf("tdps->intervals=%d, intvRadius=%d\n", tdps->intervals, intvRadius);
	
	int i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	int dataSeriesLength = r1*r2;
	//	printf ("%d %d\n", r1, r2);

	unsigned char* leadNum;

	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (double*)malloc(sizeof(double)*dataSeriesLength);

	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);

	unsigned char preBytes[8];
	unsigned char curBytes[8];

	memset(preBytes, 0, 8);

	int curByteIndex = 0;
	int reqLength, reqBytesLength, resiBitsLength, resiBits, reqExpo, reqMantLength; 
	unsigned char leadingNum;	
	double medianValue, exactData, predValue, realPrecision;
	int type_;
	double pred1D, pred2D;
	int ii, jj, II = 0, JJ = 0, updateReqLength = 1;

	int blockSize = computeBlockEdgeSize_2D(tdps->segment_size);
	int R1 = 1+(r1-1)/blockSize;
	int R2 = 1+(r2-1)/blockSize;		
	double* pwrErrBound = extractRealPrecision_2D_double(R1, R2, blockSize, tdps);

	realPrecision = pwrErrBound[0];	
	computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
	reqBytesLength = reqLength/8;
	resiBitsLength = reqLength%8;

	/* Process Row-0, data 0 */

	// compute resiBits
	resiBits = 0;
	if (resiBitsLength != 0) {
		int kMod8 = k % 8;
		int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
		if (rightMovSteps > 0) {
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
		} else if (rightMovSteps < 0) {
			int code1 = getLeftMovingCode(kMod8);
			int code2 = getRightMovingCode(kMod8, resiBitsLength);
			int leftMovSteps = -rightMovSteps;
			rightMovSteps = 8 - leftMovSteps;
			resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
			p++;
			resiBits = resiBits
					| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
		} else // rightMovSteps == 0
		{
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code);
			p++;
		}
		k += resiBitsLength;
	}

	// recover the exact data
	memset(curBytes, 0, 8);
	leadingNum = leadNum[l++];
	memcpy(curBytes, preBytes, leadingNum);
	for (j = leadingNum; j < reqBytesLength; j++)
		curBytes[j] = tdps->exactMidBytes[curByteIndex++];
	if (resiBitsLength != 0) {
		unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
		curBytes[reqBytesLength] = resiByte;
	}

	exactData = bytesToDouble(curBytes);
	(*data)[0] = exactData + medianValue;
	memcpy(preBytes,curBytes,8);

	/* Process Row-0, data 1 */
	type_ = type[1]; 
	if (type_ != 0)
	{
		pred1D = (*data)[0];		
		(*data)[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else
	{
		// compute resiBits
		resiBits = 0;
		if (resiBitsLength != 0) {
			int kMod8 = k % 8;
			int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
			if (rightMovSteps > 0) {
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
			} else if (rightMovSteps < 0) {
				int code1 = getLeftMovingCode(kMod8);
				int code2 = getRightMovingCode(kMod8, resiBitsLength);
				int leftMovSteps = -rightMovSteps;
				rightMovSteps = 8 - leftMovSteps;
				resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
				p++;
				resiBits = resiBits
						| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
			} else // rightMovSteps == 0
			{
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code);
				p++;
			}
			k += resiBitsLength;
		}

		// recover the exact data
		memset(curBytes, 0, 8);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}
		
		exactData = bytesToDouble(curBytes);
		(*data)[1] = exactData + medianValue;
		memcpy(preBytes,curBytes,8);
	}

	/* Process Row-0, data 2 --> data r2-1 */
	for (jj = 2; jj < r2; jj++)
	{
		if(jj%blockSize==0)
		{
			II = 0;
			JJ++;
			realPrecision = pwrErrBound[JJ];
			updateReqLength = 0;			
		}		
		
		type_ = type[jj];
		if (type_ != 0)
		{			
			pred1D = 2*(*data)[jj-1] - (*data)[jj-2];
			(*data)[jj] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			if(updateReqLength==0)
			{
				computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;				
				updateReqLength = 1;
			}			
			// compute resiBits
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToDouble(curBytes);
			(*data)[jj] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
		}
	}

	int index;
	/* Process Row-1 --> Row-r1-1 */
	for (ii = 1; ii < r1; ii++)
	{
		/* Process row-ii data 0 */
		if(ii%blockSize==0)
			II++;
		JJ = 0;
		realPrecision = pwrErrBound[II*R2+JJ];				
		updateReqLength = 0;

		index = ii*r2;

		type_ = type[index];
		if (type_ != 0)
		{
			pred1D = (*data)[index-r2];
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			if(updateReqLength==0)
			{
				computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;				
				updateReqLength = 1;
			}
			
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToDouble(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
		}

		/* Process row-ii data 1 --> r2-1*/
		for (jj = 1; jj < r2; jj++)
		{
			index = ii*r2+jj;

			if(jj%blockSize==0)
				JJ++;
			realPrecision = pwrErrBound[II*R2+JJ];			
			updateReqLength = 0;

			type_ = type[index];
			if (type_ != 0)
			{
				pred2D = (*data)[index-1] + (*data)[index-r2] - (*data)[index-r2-1];
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				if(updateReqLength==0)
				{
					computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
					reqBytesLength = reqLength/8;
					resiBitsLength = reqLength%8;				
					updateReqLength = 1;
				}						
				
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 8);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToDouble(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,8);
			}
		}
	}

	free(pwrErrBound);
	free(leadNum);
	free(type);
	return;
}

double* extractRealPrecision_3D_double(int R1, int R2, int R3, int blockSize, TightDataPointStorageD* tdps)
{
	int i,j,k=0, IR, JR, p = 0;
	int R23 = R2*R3;
	unsigned char* bytes = tdps->pwrErrBoundBytes;
	unsigned char tmpBytes[4] = {0};
	double* result = (double*)malloc(sizeof(double)*R1*R2*R3);
	for(i=0;i<R1;i++)
	{
		IR = i*R23;
		for(j=0;j<R2;j++)
		{
			JR = j*R3;
			for(k=0;k<R3;k++)
			{
				tmpBytes[0] = bytes[p++];
				tmpBytes[1] = bytes[p++];
				result[IR+JR+k]=bytesToDouble(tmpBytes);				
			}
		}
	}
	return result;
}

void decompressDataSeries_double_3D_pwr(double** data, int r1, int r2, int r3, TightDataPointStorageD* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	int i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	int dataSeriesLength = r1*r2*r3;
	int r23 = r2*r3;
//	printf ("%d %d %d\n", r1, r2, r3);

	unsigned char* leadNum;

	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (double*)malloc(sizeof(double)*dataSeriesLength);

	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);
	//for(i=0;i<dataSeriesLength;i++)
	//	printf("%u\n", type[i]);

	unsigned char preBytes[8];
	unsigned char curBytes[8];

	memset(preBytes, 0, 8);

	int curByteIndex = 0;
	int reqLength, reqBytesLength, resiBitsLength, resiBits, reqExpo, reqMantLength; 
	unsigned char leadingNum;
	double medianValue, exactData, predValue, realPrecision;
	int type_;
	double pred1D, pred2D, pred3D;
	int ii, jj, kk, II = 0, JJ = 0, KK = 0, updateReqLength = 1;

	int blockSize = computeBlockEdgeSize_3D(tdps->segment_size);
	int R1 = 1+(r1-1)/blockSize;
	int R2 = 1+(r2-1)/blockSize;		
	int R3 = 1+(r3-1)/blockSize;
	int R23 = R2*R3;
	double* pwrErrBound = extractRealPrecision_3D_double(R1, R2, R3, blockSize, tdps);

	realPrecision = pwrErrBound[0];	
	computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
	reqBytesLength = reqLength/8;
	resiBitsLength = reqLength%8;

	///////////////////////////	Process layer-0 ///////////////////////////
	/* Process Row-0 data 0*/
	// compute resiBits
	resiBits = 0;
	if (resiBitsLength != 0) {
		int kMod8 = k % 8;
		int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
		if (rightMovSteps > 0) {
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
		} else if (rightMovSteps < 0) {
			int code1 = getLeftMovingCode(kMod8);
			int code2 = getRightMovingCode(kMod8, resiBitsLength);
			int leftMovSteps = -rightMovSteps;
			rightMovSteps = 8 - leftMovSteps;
			resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
			p++;
			resiBits = resiBits
					| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
		} else // rightMovSteps == 0
		{
			int code = getRightMovingCode(kMod8, resiBitsLength);
			resiBits = (tdps->residualMidBits[p] & code);
			p++;
		}
		k += resiBitsLength;
	}

	// recover the exact data
	memset(curBytes, 0, 8);
	leadingNum = leadNum[l++];
	memcpy(curBytes, preBytes, leadingNum);
	for (j = leadingNum; j < reqBytesLength; j++)
		curBytes[j] = tdps->exactMidBytes[curByteIndex++];
	if (resiBitsLength != 0) {
		unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
		curBytes[reqBytesLength] = resiByte;
	}

	exactData = bytesToDouble(curBytes);
	(*data)[0] = exactData + medianValue;
	memcpy(preBytes,curBytes,8);

	/* Process Row-0, data 1 */
	pred1D = (*data)[0];

	type_ = type[1];
	if (type_ != 0)
	{
		(*data)[1] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
	}
	else
	{
		// compute resiBits
		resiBits = 0;
		if (resiBitsLength != 0) {
			int kMod8 = k % 8;
			int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
			if (rightMovSteps > 0) {
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
			} else if (rightMovSteps < 0) {
				int code1 = getLeftMovingCode(kMod8);
				int code2 = getRightMovingCode(kMod8, resiBitsLength);
				int leftMovSteps = -rightMovSteps;
				rightMovSteps = 8 - leftMovSteps;
				resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
				p++;
				resiBits = resiBits
						| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
			} else // rightMovSteps == 0
			{
				int code = getRightMovingCode(kMod8, resiBitsLength);
				resiBits = (tdps->residualMidBits[p] & code);
				p++;
			}
			k += resiBitsLength;
		}

		// recover the exact data
		memset(curBytes, 0, 8);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}

		exactData = bytesToDouble(curBytes);
		(*data)[1] = exactData + medianValue;
		memcpy(preBytes,curBytes,8);
	}

	/* Process Row-0, data 2 --> data r3-1 */
	for (jj = 2; jj < r3; jj++)
	{
		if(jj%blockSize==0)
		{
			KK = 0;//dimension 1 (top)
			II = 0;//dimension 2 (mid)
			JJ++;
			realPrecision = pwrErrBound[JJ];
			updateReqLength = 0;			
		}		
		type_ = type[jj];
		if (type_ != 0)
		{
			pred1D = 2*(*data)[jj-1] - (*data)[jj-2];
			(*data)[jj] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			if(updateReqLength==0)
			{
				computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;				
				updateReqLength = 1;
			}

			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToDouble(curBytes);
			(*data)[jj] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
		}
	}

	int index;
	/* Process Row-1 --> Row-r2-1 */
	for (ii = 1; ii < r2; ii++)
	{
		/* Process row-ii data 0 */		
		if(ii%blockSize==0)
			II++;		
		JJ = 0;
		realPrecision = pwrErrBound[II*R3+JJ];
		updateReqLength = 0;		

		index = ii*r3;
		
		type_ = type[index];
		if (type_ != 0)
		{
			pred1D = (*data)[index-r3];			
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			if(updateReqLength==0)
			{
				computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;				
				updateReqLength = 1;
			}
			
			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToDouble(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
		}

		/* Process row-ii data 1 --> r3-1*/
		for (jj = 1; jj < r3; jj++)
		{
			index = ii*r3+jj;

			if(jj%blockSize==0)
				JJ++;
			realPrecision = pwrErrBound[II*R3+JJ];			
			updateReqLength = 0;			
			
			type_ = type[index];
			if (type_ != 0)
			{
				pred2D = (*data)[index-1] + (*data)[index-r3] - (*data)[index-r3-1];				
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				if(updateReqLength==0)
				{
					computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
					reqBytesLength = reqLength/8;
					resiBitsLength = reqLength%8;				
					updateReqLength = 1;
				}

				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 8);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToDouble(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,8);
			}
		}
	}

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (kk = 1; kk < r1; kk++)
	{
		/* Process Row-0 data 0*/
		index = kk*r23;		
		if(kk%blockSize==0)
			KK++;
		II = 0;
		JJ = 0;

		realPrecision = pwrErrBound[KK*R23];			
		updateReqLength = 0;			

		type_ = type[index];
		if (type_ != 0)
		{
			pred1D = (*data)[index-r23];			
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
		}
		else
		{
			// compute resiBits
			if(updateReqLength==0)
			{
				computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;				
				updateReqLength = 1;
			}

			resiBits = 0;
			if (resiBitsLength != 0) {
				int kMod8 = k % 8;
				int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
				if (rightMovSteps > 0) {
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
				} else if (rightMovSteps < 0) {
					int code1 = getLeftMovingCode(kMod8);
					int code2 = getRightMovingCode(kMod8, resiBitsLength);
					int leftMovSteps = -rightMovSteps;
					rightMovSteps = 8 - leftMovSteps;
					resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
					p++;
					resiBits = resiBits
							| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
				} else // rightMovSteps == 0
				{
					int code = getRightMovingCode(kMod8, resiBitsLength);
					resiBits = (tdps->residualMidBits[p] & code);
					p++;
				}
				k += resiBitsLength;
			}

			// recover the exact data
			memset(curBytes, 0, 8);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToDouble(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,8);
		}

		/* Process Row-0 data 1 --> data r3-1 */
		for (jj = 1; jj < r3; jj++)
		{
			index = kk*r23+jj;

			if(index==46788)
				printf("index=%d\n", index);
			if(jj%blockSize==0)
				JJ++;

			realPrecision = pwrErrBound[KK*R23+JJ];			
			updateReqLength = 0;			
			
			type_ = type[index];
			if (type_ != 0)
			{
				pred2D = (*data)[index-1] + (*data)[index-r23] - (*data)[index-r23-1];			
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				if(updateReqLength==0)
				{
					computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
					reqBytesLength = reqLength/8;
					resiBitsLength = reqLength%8;				
					updateReqLength = 1;
				}
			
				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 8);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToDouble(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,8);
			}
		}

		/* Process Row-1 --> Row-r2-1 */
		for (ii = 1; ii < r2; ii++)
		{
			/* Process Row-i data 0 */
			index = kk*r23 + ii*r3;
			
			if(ii%blockSize==0)
				II++;
			JJ = 0;
			
			realPrecision = pwrErrBound[KK*R23+II*R3];			
			updateReqLength = 0;						

			type_ = type[index];
			if (type_ != 0)
			{
				pred2D = (*data)[index-r3] + (*data)[index-r23] - (*data)[index-r23-r3];				
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
			}
			else
			{
				// compute resiBits
				if(updateReqLength==0)
				{
					computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
					reqBytesLength = reqLength/8;
					resiBitsLength = reqLength%8;				
					updateReqLength = 1;
				}

				resiBits = 0;
				if (resiBitsLength != 0) {
					int kMod8 = k % 8;
					int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
					if (rightMovSteps > 0) {
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
					} else if (rightMovSteps < 0) {
						int code1 = getLeftMovingCode(kMod8);
						int code2 = getRightMovingCode(kMod8, resiBitsLength);
						int leftMovSteps = -rightMovSteps;
						rightMovSteps = 8 - leftMovSteps;
						resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
						p++;
						resiBits = resiBits
								| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
					} else // rightMovSteps == 0
					{
						int code = getRightMovingCode(kMod8, resiBitsLength);
						resiBits = (tdps->residualMidBits[p] & code);
						p++;
					}
					k += resiBitsLength;
				}

				// recover the exact data
				memset(curBytes, 0, 8);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToDouble(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,8);
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (jj = 1; jj < r3; jj++)
			{
				index = kk*r23 + ii*r3 + jj;
				if(jj%blockSize==0)
					JJ++;

				realPrecision = pwrErrBound[KK*R23+II*R3+JJ];			
				updateReqLength = 0;				

				type_ = type[index];
				if (type_ != 0)
				{
					pred3D = (*data)[index-1] + (*data)[index-r3] + (*data)[index-r23]
					- (*data)[index-r3-1] - (*data)[index-r23-r3] - (*data)[index-r23-1] + (*data)[index-r23-r3-1];					
					(*data)[index] = pred3D + 2 * (type_ - intvRadius) * realPrecision;
				}
				else
				{
					// compute resiBits
					if(updateReqLength==0)
					{
						computeReqLength_double(realPrecision, tdps->radExpo, &reqLength, &medianValue);
						reqBytesLength = reqLength/8;
						resiBitsLength = reqLength%8;				
						updateReqLength = 1;
					}
				
					resiBits = 0;
					if (resiBitsLength != 0) {
						int kMod8 = k % 8;
						int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
						if (rightMovSteps > 0) {
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
						} else if (rightMovSteps < 0) {
							int code1 = getLeftMovingCode(kMod8);
							int code2 = getRightMovingCode(kMod8, resiBitsLength);
							int leftMovSteps = -rightMovSteps;
							rightMovSteps = 8 - leftMovSteps;
							resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
							p++;
							resiBits = resiBits
									| ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
						} else // rightMovSteps == 0
						{
							int code = getRightMovingCode(kMod8, resiBitsLength);
							resiBits = (tdps->residualMidBits[p] & code);
							p++;
						}
						k += resiBitsLength;
					}

					// recover the exact data
					memset(curBytes, 0, 8);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToDouble(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,8);
				}
			}
		}
	}

	free(pwrErrBound);
	free(leadNum);
	free(type);
	return;
}

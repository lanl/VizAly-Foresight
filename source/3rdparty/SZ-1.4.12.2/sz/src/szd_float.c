/**
 *  @file szd_float.c
 *  @author Sheng Di and Dingwen Tao
 *  @date Aug, 2016
 *  @brief 
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include "szd_float.h"
#include "TightDataPointStorageF.h"
#include "sz.h"
#include "Huffman.h"
#include "szd_float_pwr.h"
//#include "rw.h"

/**
 * 
 * 
 * @return status SUCCESSFUL (SZ_SCES) or not (other error codes) f
 * */
int SZ_decompress_args_float(float** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize)
{
	int status = SZ_SCES;
	size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
	
	//unsigned char* tmpBytes;
	size_t targetUncompressSize = dataLength <<2; //i.e., *4
	//tmpSize must be "much" smaller than dataLength
	size_t i, tmpSize = 8+MetaDataByteLength+SZ_SIZE_TYPE;
	unsigned char* szTmpBytes;	
	
	if(cmpSize!=8+4+MetaDataByteLength && cmpSize!=8+8+MetaDataByteLength) //4,8 means two posibilities of SZ_SIZE_TYPE
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
			tmpSize = zlib_uncompress5(cmpBytes, (unsigned long)cmpSize, &szTmpBytes, (unsigned long)targetUncompressSize+4+MetaDataByteLength+SZ_SIZE_TYPE);//		(unsigned long)targetUncompressSize+8: consider the total length under lossless compression mode is actually 3+4+1+targetUncompressSize
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
			memcpy(*newData, szTmpBytes+4+MetaDataByteLength+SZ_SIZE_TYPE, dataLength*floatSize);
		}
		else
		{
			unsigned char* p = szTmpBytes+4+MetaDataByteLength+SZ_SIZE_TYPE;
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
	if(szMode!=SZ_BEST_SPEED && cmpSize!=8+MetaDataByteLength+SZ_SIZE_TYPE)
		free(szTmpBytes);
	SZ_ReleaseHuffman();	
	return status;
}

void decompressDataSeries_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	size_t i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
								// in resiMidBits, p is to track the
								// byte_index of resiMidBits, l is for
								// leadNum
	unsigned char* leadNum;
	double interval = tdps->realPrecision*2;
	
	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (float*)malloc(sizeof(float)*dataSeriesLength);

	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);

	//sdi:Debug
	//writeUShortData(type, dataSeriesLength, "decompressStateBytes.sb");

	unsigned char preBytes[4];
	unsigned char curBytes[4];
	
	memset(preBytes, 0, 4);

	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits; 
	unsigned char leadingNum;	
	float medianValue, exactData, predValue;
	
	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;
	
	int type_;
	for (i = 0; i < dataSeriesLength; i++) {
		type_ = type[i];
		switch (type_) {
		case 0:
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
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}
			
			exactData = bytesToFloat(curBytes);
			(*data)[i] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
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

void decompressDataSeries_float_2D(float** data, size_t r1, size_t r2, TightDataPointStorageF* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	//printf("tdps->intervals=%d, intvRadius=%d\n", tdps->intervals, intvRadius);
	
	size_t i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	size_t dataSeriesLength = r1*r2;
	//	printf ("%d %d\n", r1, r2);

	unsigned char* leadNum;
	double realPrecision = tdps->realPrecision;

	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (float*)malloc(sizeof(float)*dataSeriesLength);

	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);

	unsigned char preBytes[4];
	unsigned char curBytes[4];

	memset(preBytes, 0, 4);

	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits; 
	unsigned char leadingNum;	
	float medianValue, exactData, predValue;
	int type_;

	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;
	
	float pred1D, pred2D;
	size_t ii, jj;

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
	memset(curBytes, 0, 4);
	leadingNum = leadNum[l++];
	memcpy(curBytes, preBytes, leadingNum);
	for (j = leadingNum; j < reqBytesLength; j++)
		curBytes[j] = tdps->exactMidBytes[curByteIndex++];
	if (resiBitsLength != 0) {
		unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
		curBytes[reqBytesLength] = resiByte;
	}

	exactData = bytesToFloat(curBytes);
	(*data)[0] = exactData + medianValue;
	memcpy(preBytes,curBytes,4);

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
		memset(curBytes, 0, 4);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}

		exactData = bytesToFloat(curBytes);
		(*data)[1] = exactData + medianValue;
		memcpy(preBytes,curBytes,4);
	}

	/* Process Row-0, data 2 --> data r2-1 */
	for (jj = 2; jj < r2; jj++)
	{
		type_ = type[jj];
		if (type_ != 0)
		{
			pred1D = 2*(*data)[jj-1] - (*data)[jj-2];				
			(*data)[jj] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
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
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[jj] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}
	}

	size_t index;
	/* Process Row-1 --> Row-r1-1 */
	for (ii = 1; ii < r1; ii++)
	{
		/* Process row-ii data 0 */
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
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}

		/* Process row-ii data 1 --> r2-1*/
		for (jj = 1; jj < r2; jj++)
		{
			index = ii*r2+jj;
			pred2D = (*data)[index-1] + (*data)[index-r2] - (*data)[index-r2-1];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
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
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}
		}
	}

	free(leadNum);
	free(type);
	return;
}

void decompressDataSeries_float_3D(float** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageF* tdps) 
{
	updateQuantizationInfo(tdps->intervals);
	size_t i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	size_t dataSeriesLength = r1*r2*r3;
	size_t r23 = r2*r3;
//	printf ("%d %d %d\n", r1, r2, r3);
	unsigned char* leadNum;
	double realPrecision = tdps->realPrecision;

	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (float*)malloc(sizeof(float)*dataSeriesLength);
	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);

	unsigned char preBytes[4];
	unsigned char curBytes[4];

	memset(preBytes, 0, 4);
	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits;
	unsigned char leadingNum;
	float medianValue, exactData, predValue;
	int type_;

	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;
	
	float pred1D, pred2D, pred3D;
	size_t ii, jj, kk;

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
	memset(curBytes, 0, 4);
	leadingNum = leadNum[l++];
	memcpy(curBytes, preBytes, leadingNum);
	for (j = leadingNum; j < reqBytesLength; j++)
		curBytes[j] = tdps->exactMidBytes[curByteIndex++];
	if (resiBitsLength != 0) {
		unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
		curBytes[reqBytesLength] = resiByte;
	}
	exactData = bytesToFloat(curBytes);
	(*data)[0] = exactData + medianValue;
	memcpy(preBytes,curBytes,4);

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
		memset(curBytes, 0, 4);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}

		exactData = bytesToFloat(curBytes);
		(*data)[1] = exactData + medianValue;
		memcpy(preBytes,curBytes,4);
	}
	/* Process Row-0, data 2 --> data r3-1 */
	for (jj = 2; jj < r3; jj++)
	{
		pred1D = 2*(*data)[jj-1] - (*data)[jj-2];

		type_ = type[jj];
		if (type_ != 0)
		{
			(*data)[jj] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
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
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[jj] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}
	}

	size_t index;
	/* Process Row-1 --> Row-r2-1 */
	for (ii = 1; ii < r2; ii++)
	{
		/* Process row-ii data 0 */
		index = ii*r3;
		pred1D = (*data)[index-r3];

		type_ = type[index];
		if (type_ != 0)
		{
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
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
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}

		/* Process row-ii data 1 --> r3-1*/
		for (jj = 1; jj < r3; jj++)
		{
			index = ii*r3+jj;
			pred2D = (*data)[index-1] + (*data)[index-r3] - (*data)[index-r3-1];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
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
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}
		}
	}

	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (kk = 1; kk < r1; kk++)
	{
		/* Process Row-0 data 0*/
		index = kk*r23;
		pred1D = (*data)[index-r23];

		type_ = type[index];
		if (type_ != 0)
		{
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
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
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}

		/* Process Row-0 data 1 --> data r3-1 */
		for (jj = 1; jj < r3; jj++)
		{
			index = kk*r23+jj;
			pred2D = (*data)[index-1] + (*data)[index-r23] - (*data)[index-r23-1];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
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
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}
		}

		/* Process Row-1 --> Row-r2-1 */
		for (ii = 1; ii < r2; ii++)
		{
			/* Process Row-i data 0 */
			index = kk*r23 + ii*r3;
			pred2D = (*data)[index-r3] + (*data)[index-r23] - (*data)[index-r23-r3];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
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
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (jj = 1; jj < r3; jj++)
			{
				index = kk*r23 + ii*r3 + jj;
				pred3D = (*data)[index-1] + (*data)[index-r3] + (*data)[index-r23]
					- (*data)[index-r3-1] - (*data)[index-r23-r3] - (*data)[index-r23-1] + (*data)[index-r23-r3-1];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred3D + 2 * (type_ - intvRadius) * realPrecision;
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
					memset(curBytes, 0, 4);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToFloat(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,4);
				}
			}
		}
	}

	free(leadNum);
	free(type);
	return;
}


void decompressDataSeries_float_4D(float** data, size_t r1, size_t r2, size_t r3, size_t r4, TightDataPointStorageF* tdps)
{
	updateQuantizationInfo(tdps->intervals);
	size_t i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
	// in resiMidBits, p is to track the
	// byte_index of resiMidBits, l is for
	// leadNum
	size_t dataSeriesLength = r1*r2*r3*r4;
	size_t r234 = r2*r3*r4;
	size_t r34 = r3*r4;
//	printf ("%d %d %d %d\n", r1, r2, r3, r4);
	unsigned char* leadNum;
	double realPrecision = tdps->realPrecision;

	convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

	*data = (float*)malloc(sizeof(float)*dataSeriesLength);
	int* type = (int*)malloc(dataSeriesLength*sizeof(int));
	//convertByteArray2IntArray_fast_3b(dataSeriesLength, tdps->typeArray, tdps->typeArray_size, &type);
	//reconstruct_HuffTree_and_Decode_16states(tdps->typeArray, dataSeriesLength, &type);
	//memcpy(type, tdps->typeArray, dataSeriesLength*sizeof(unsigned short));
	//type = tdps->typeArray;
	decode_withTree(tdps->typeArray, dataSeriesLength, type);

	unsigned char preBytes[4];
	unsigned char curBytes[4];

	memset(preBytes, 0, 4);
	size_t curByteIndex = 0;
	int reqBytesLength, resiBitsLength, resiBits;
	unsigned char leadingNum;
	float medianValue, exactData, predValue;
	int type_;

	reqBytesLength = tdps->reqLength/8;
	resiBitsLength = tdps->reqLength%8;
	medianValue = tdps->medianValue;

	float pred1D, pred2D, pred3D;
	size_t ii, jj, kk, ll;
	size_t index;

	for (ll = 0; ll < r1; ll++)
	{

		///////////////////////////	Process layer-0 ///////////////////////////
		/* Process Row-0 data 0*/
		index = ll*r234;

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
		memset(curBytes, 0, 4);
		leadingNum = leadNum[l++];
		memcpy(curBytes, preBytes, leadingNum);
		for (j = leadingNum; j < reqBytesLength; j++)
			curBytes[j] = tdps->exactMidBytes[curByteIndex++];
		if (resiBitsLength != 0) {
			unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
			curBytes[reqBytesLength] = resiByte;
		}
		exactData = bytesToFloat(curBytes);
		(*data)[index] = exactData + medianValue;
		memcpy(preBytes,curBytes,4);

		/* Process Row-0, data 1 */
		index = ll*r234+1;

		pred1D = (*data)[index-1];

		type_ = type[index];
		if (type_ != 0)
		{
			(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
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
			memset(curBytes, 0, 4);
			leadingNum = leadNum[l++];
			memcpy(curBytes, preBytes, leadingNum);
			for (j = leadingNum; j < reqBytesLength; j++)
				curBytes[j] = tdps->exactMidBytes[curByteIndex++];
			if (resiBitsLength != 0) {
				unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
				curBytes[reqBytesLength] = resiByte;
			}

			exactData = bytesToFloat(curBytes);
			(*data)[index] = exactData + medianValue;
			memcpy(preBytes,curBytes,4);
		}

		/* Process Row-0, data 2 --> data r4-1 */
		for (jj = 2; jj < r4; jj++)
		{
			index = ll*r234+jj;

			pred1D = 2*(*data)[index-1] - (*data)[index-2];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
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
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}
		}

		/* Process Row-1 --> Row-r3-1 */
		for (ii = 1; ii < r3; ii++)
		{
			/* Process row-ii data 0 */
			index = ll*r234+ii*r4;

			pred1D = (*data)[index-r4];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
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
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}

			/* Process row-ii data 1 --> r4-1*/
			for (jj = 1; jj < r4; jj++)
			{
				index = ll*r234+ii*r4+jj;

				pred2D = (*data)[index-1] + (*data)[index-r4] - (*data)[index-r4-1];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
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
					memset(curBytes, 0, 4);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToFloat(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,4);
				}
			}
		}

		///////////////////////////	Process layer-1 --> layer-r2-1 ///////////////////////////

		for (kk = 1; kk < r2; kk++)
		{
			/* Process Row-0 data 0*/
			index = ll*r234+kk*r34;

			pred1D = (*data)[index-r34];

			type_ = type[index];
			if (type_ != 0)
			{
				(*data)[index] = pred1D + 2 * (type_ - intvRadius) * realPrecision;
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
				memset(curBytes, 0, 4);
				leadingNum = leadNum[l++];
				memcpy(curBytes, preBytes, leadingNum);
				for (j = leadingNum; j < reqBytesLength; j++)
					curBytes[j] = tdps->exactMidBytes[curByteIndex++];
				if (resiBitsLength != 0) {
					unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
					curBytes[reqBytesLength] = resiByte;
				}

				exactData = bytesToFloat(curBytes);
				(*data)[index] = exactData + medianValue;
				memcpy(preBytes,curBytes,4);
			}

			/* Process Row-0 data 1 --> data r4-1 */
			for (jj = 1; jj < r4; jj++)
			{
				index = ll*r234+kk*r34+jj;

				pred2D = (*data)[index-1] + (*data)[index-r34] - (*data)[index-r34-1];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
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
					memset(curBytes, 0, 4);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToFloat(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,4);
				}
			}

			/* Process Row-1 --> Row-r3-1 */
			for (ii = 1; ii < r3; ii++)
			{
				/* Process Row-i data 0 */
				index = ll*r234+kk*r34+ii*r4;

				pred2D = (*data)[index-r4] + (*data)[index-r34] - (*data)[index-r34-r4];

				type_ = type[index];
				if (type_ != 0)
				{
					(*data)[index] = pred2D + 2 * (type_ - intvRadius) * realPrecision;
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
					memset(curBytes, 0, 4);
					leadingNum = leadNum[l++];
					memcpy(curBytes, preBytes, leadingNum);
					for (j = leadingNum; j < reqBytesLength; j++)
						curBytes[j] = tdps->exactMidBytes[curByteIndex++];
					if (resiBitsLength != 0) {
						unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
						curBytes[reqBytesLength] = resiByte;
					}

					exactData = bytesToFloat(curBytes);
					(*data)[index] = exactData + medianValue;
					memcpy(preBytes,curBytes,4);
				}

				/* Process Row-i data 1 --> data r4-1 */
				for (jj = 1; jj < r4; jj++)
				{
					index = ll*r234+kk*r34+ii*r4+jj;

					pred3D = (*data)[index-1] + (*data)[index-r4] + (*data)[index-r34]
							- (*data)[index-r4-1] - (*data)[index-r34-r4] - (*data)[index-r34-1] + (*data)[index-r34-r4-1];


					type_ = type[index];
					if (type_ != 0)
					{
						(*data)[index] = pred3D + 2 * (type_ - intvRadius) * realPrecision;
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
						memset(curBytes, 0, 4);
						leadingNum = leadNum[l++];
						memcpy(curBytes, preBytes, leadingNum);
						for (j = leadingNum; j < reqBytesLength; j++)
							curBytes[j] = tdps->exactMidBytes[curByteIndex++];
						if (resiBitsLength != 0) {
							unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
							curBytes[reqBytesLength] = resiByte;
						}

						exactData = bytesToFloat(curBytes);
						(*data)[index] = exactData + medianValue;
						memcpy(preBytes,curBytes,4);
					}
				}
			}

		}
	}

	free(leadNum);
	free(type);
	return;
}

void getSnapshotData_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps, int errBoundMode)
{	
	SZ_Reset();
	size_t i;

	if (tdps->allSameData) {
		float value = bytesToFloat(tdps->exactMidBytes);
		*data = (float*)malloc(sizeof(float)*dataSeriesLength);
		for (i = 0; i < dataSeriesLength; i++)
			(*data)[i] = value;
	} else {
		if (tdps->rtypeArray == NULL) {
			if(errBoundMode < PW_REL)
				decompressDataSeries_float_1D(data, dataSeriesLength, tdps);
			else 
			{
				//decompressDataSeries_float_1D_pwr(data, dataSeriesLength, tdps);
				decompressDataSeries_float_1D_pwrgroup(data, dataSeriesLength, tdps);
			}
			return;
		} else {
			*data = (float*)malloc(sizeof(float)*dataSeriesLength);
			// insert the reserved values
			//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
			//		dataSeriesLength, rtypeArray);
			int* rtypes;
			int validLength = computeBitNumRequired(dataSeriesLength);
			decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
			size_t count = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 1)
					(*data)[i] = tdps->reservedValue;
				else
					count++;
			}
			// get the decompressed data
			float* decmpData;
			if(errBoundMode < PW_REL)
				decompressDataSeries_float_1D(&decmpData, dataSeriesLength, tdps);
			else 
				decompressDataSeries_float_1D_pwr(&decmpData, dataSeriesLength, tdps);
			// insert the decompressed data
			size_t k = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 0) {
					(*data)[i] = decmpData[k++];
				}
			}
			free(decmpData);
			free(rtypes);
		}
	}
}

void getSnapshotData_float_2D(float** data, size_t r1, size_t r2, TightDataPointStorageF* tdps, int errBoundMode) 
{
	SZ_Reset();
	size_t i;
	size_t dataSeriesLength = r1*r2;
	if (tdps->allSameData) {
		float value = bytesToFloat(tdps->exactMidBytes);
		*data = (float*)malloc(sizeof(float)*dataSeriesLength);
		for (i = 0; i < dataSeriesLength; i++)
			(*data)[i] = value;
	} else {
		if (tdps->rtypeArray == NULL) {
			if(errBoundMode < PW_REL)
				decompressDataSeries_float_2D(data, r1, r2, tdps);
			else 
				decompressDataSeries_float_2D_pwr(data, r1, r2, tdps);
			return;
		} else {
			*data = (float*)malloc(sizeof(float)*dataSeriesLength);
			// insert the reserved values
			//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
			//		dataSeriesLength, rtypeArray);
			int* rtypes;
			int validLength = computeBitNumRequired(dataSeriesLength);
			decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
			size_t count = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 1)
					(*data)[i] = tdps->reservedValue;
				else
					count++;
			}
			// get the decompressed data
			float* decmpData;
			if(errBoundMode < PW_REL)
				decompressDataSeries_float_2D(&decmpData, r1, r2, tdps);
			else 
				decompressDataSeries_float_2D_pwr(&decmpData, r1, r2, tdps);
			// insert the decompressed data
			size_t k = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 0) {
					(*data)[i] = decmpData[k++];
				}
			}
			free(decmpData);
			free(rtypes);
		}
	}
}

void getSnapshotData_float_3D(float** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageF* tdps, int errBoundMode)
{
	SZ_Reset();
	size_t i;
	size_t dataSeriesLength = r1*r2*r3;
	if (tdps->allSameData) {
		float value = bytesToFloat(tdps->exactMidBytes);
		*data = (float*)malloc(sizeof(float)*dataSeriesLength);
		for (i = 0; i < dataSeriesLength; i++)
			(*data)[i] = value;
	} else {
		if (tdps->rtypeArray == NULL) {
			if(errBoundMode < PW_REL)
				decompressDataSeries_float_3D(data, r1, r2, r3, tdps);
			else 
				decompressDataSeries_float_3D_pwr(data, r1, r2, r3, tdps);
			return;
		} else {
			*data = (float*)malloc(sizeof(float)*dataSeriesLength);
			// insert the reserved values
			//int[] rtypes = TypeManager.convertByteArray2IntArray_fast_1b(
			//		dataSeriesLength, rtypeArray);
			int* rtypes;
			int validLength = computeBitNumRequired(dataSeriesLength);
			decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
			size_t count = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 1)
					(*data)[i] = tdps->reservedValue;
				else
					count++;
			}
			// get the decompressed data
			float* decmpData;
			if(errBoundMode < PW_REL)
				decompressDataSeries_float_3D(&decmpData, r1, r2, r3, tdps);
			else 
				decompressDataSeries_float_3D_pwr(&decmpData, r1, r2, r3, tdps);
			// insert the decompressed data
			size_t k = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 0) {
					(*data)[i] = decmpData[k++];
				}
			}
			free(decmpData);
			free(rtypes);
		}
	}
}

void getSnapshotData_float_4D(float** data, size_t r1, size_t r2, size_t r3, size_t r4, TightDataPointStorageF* tdps, int errBoundMode)
{
	SZ_Reset();
	size_t i;
	size_t dataSeriesLength = r1*r2*r3*r4;
	if (tdps->allSameData) {
		float value = bytesToFloat(tdps->exactMidBytes);
		*data = (float*)malloc(sizeof(float)*dataSeriesLength);
		for (i = 0; i < dataSeriesLength; i++)
			(*data)[i] = value;
	} else {
		if (tdps->rtypeArray == NULL) {
			if(errBoundMode < PW_REL)
				decompressDataSeries_float_4D(data, r1, r2, r3, r4, tdps);
			else
				decompressDataSeries_float_3D_pwr(data, r1*r2, r3, r4, tdps);
				//ToDO
				//decompressDataSeries_float_4D_pwr(data, r1, r2, r3, r4, tdps);
			return;
		} else {
			*data = (float*)malloc(sizeof(float)*dataSeriesLength);
			int* rtypes;
			int validLength = computeBitNumRequired(dataSeriesLength);
			decompressBitArraybySimpleLZ77(&rtypes, tdps->rtypeArray, tdps->rtypeArray_size, dataSeriesLength, validLength);
			size_t count = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 1)
					(*data)[i] = tdps->reservedValue;
				else
					count++;
			}
			// get the decompressed data
			float* decmpData;
			if(errBoundMode < PW_REL)
				decompressDataSeries_float_4D(&decmpData, r1, r2, r3, r4, tdps);
			else
				decompressDataSeries_float_3D_pwr(&decmpData, r1*r2, r3, r4, tdps);
				//ToDO
				//decompressDataSeries_float_4D_pwr(&decompData, r1, r2, r3, r4, tdps);
			// insert the decompressed data
			size_t k = 0;
			for (i = 0; i < dataSeriesLength; i++) {
				if (rtypes[i] == 0) {
					(*data)[i] = decmpData[k++];
				}
			}
			free(decmpData);
			free(rtypes);
		}
	}
}

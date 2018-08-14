/**
 *  @file sz.c
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
#include "sz.h"
#include "CompressElement.h"
#include "DynamicByteArray.h"
#include "DynamicIntArray.h"
#include "TightDataPointStorageD.h"
#include "TightDataPointStorageF.h"
#include "zlib.h"
#include "rw.h"
//#include "CurveFillingCompressStorage.h"

unsigned int maxRangeRadius = 32768;

int sysEndianType; //endian type of the system
int dataEndianType; //endian type of the data

char maxHeap[10];

long status;

int sol_ID;
int errorBoundMode; //ABS, REL, ABS_AND_REL, or ABS_OR_REL, or PW_REL

int gzipMode; //four options: Z_NO_COMPRESSION, or Z_BEST_SPEED, Z_BEST_COMPRESSION, Z_DEFAULT_COMPRESSION

char *sz_cfgFile;

int offset;

double absErrBound;
double relBoundRatio;
double psnr;
double pw_relBoundRatio;
int segment_size;
int pwr_type = SZ_PWR_MIN_TYPE;

int versionNumber[4] = {SZ_VER_MAJOR,SZ_VER_MINOR,SZ_VER_BUILD,SZ_VER_REVISION};

int spaceFillingCurveTransform; //default is 0, or 1 set by sz.config
int reOrgSize; //the granularity of the reganization of the original data

int intvCapacity = 0;
int intvRadius = 0;

int layers = 1;
float predThreshold = 0.98;
int sampleDistance = 10;
char optQuantMode = 0; //opt Quantization (0: fixed ; 1: optimized)

int szMode = SZ_BEST_COMPRESSION;

int SZ_SIZE_TYPE = 8;

SZ_VarSet* sz_varset = NULL;

sz_params *conf_params = NULL;

int SZ_Init(char *configFilePath)
{
	char str[512]="", str2[512]="", str3[512]="";
	sz_cfgFile = configFilePath;
	int loadFileResult = SZ_LoadConf();
	if(loadFileResult==SZ_NSCS)
		return SZ_NSCS;
	
	SZ_SIZE_TYPE = sizeof(size_t);
	return SZ_SCES;
}

void SZ_Reset()
{
    if(pool==NULL)
    {
		pool = (struct node_t*)malloc(allNodes*2*sizeof(struct node_t));
		qqq = (node*)malloc(allNodes*2*sizeof(node));
		code = (unsigned long**)malloc(stateNum*sizeof(unsigned long*));
		cout = (unsigned char *)malloc(stateNum*sizeof(unsigned char));
	}
	
	memset(pool, 0, allNodes*2*sizeof(struct node_t));
	memset(qqq, 0, allNodes*2*sizeof(node));
    memset(code, 0, stateNum*sizeof(unsigned long*));
    memset(cout, 0, stateNum*sizeof(unsigned char));
	qq = qqq - 1;
	n_nodes = 0;
    n_inode = 0;
    qend = 1;
}

int SZ_Init_Params(sz_params *params)
{
	conf_params = (sz_params*)malloc(sizeof(sz_params));   
	memcpy(conf_params, params, sizeof(sz_params));
    int x = 1;
    char *y = (char*)&x;
    int endianType = BIG_ENDIAN_SYSTEM;
    if(*y==1) endianType = LITTLE_ENDIAN_SYSTEM;

	SZ_SIZE_TYPE = sizeof(size_t);

    // set default values
    if(params->max_quant_intervals > 0) 
		maxRangeRadius = params->max_quant_intervals/2;
	else
		params->max_quant_intervals = maxRangeRadius*2;

	stateNum = maxRangeRadius*2;
	allNodes = maxRangeRadius*4;
	intvCapacity = maxRangeRadius*2;
	intvRadius = maxRangeRadius;

    dataEndianType    = endianType;
    //sysEndianType    = endianType;
    sol_ID                    = SZ;
    offset                    = 0;
    gzipMode               = Z_BEST_SPEED;
    sampleDistance = 50;
    predThreshold = 0.97;
    errorBoundMode    = REL;
    absErrBound         = 0.000001;
    relBoundRatio         = 0.001;
    szMode = SZ_BEST_COMPRESSION;

    // set values from function arguments if avail.
    // [ENV]
    if(params->dataEndianType >= 0) dataEndianType    = params->dataEndianType;
    //if(params->sysEndianType >= 0)    sysEndianType    = params->sysEndianType;
    if(params->sol_ID >= 0)  
		sol_ID = params->sol_ID;

    // [PARAMETER]
    if(sol_ID==SZ) {
        if(params->offset >= 0) offset = params->offset;

        /* gzipModes:
            Gzip_NO_COMPRESSION=0,
            Gzip_BEST_SPEED=1,
            Gzip_BEST_COMPRESSION=9,
            Gzip_DEFAULT_COMPRESSION=-1 */
        if(params->gzipMode >= -1) gzipMode = params->gzipMode;

		if(params->szMode >= 0) szMode = params->szMode;
        //if(params->maxSegmentNum >= 0) maxSegmentNum = params->maxSegmentNum;
        //if(params->spaceFillingCurveTransform >= 0) spaceFillingCurveTransform = params->spaceFillingCurveTransform;
        //if(params->reOrgSize >= 0) reOrgSize = params->reOrgSize;

        /* errBoundModes:
            ABS = 0
            REL = 1
            ABS_AND_REL =2
            ABS_OR_REL = 3 */
        if( params->errorBoundMode >= 0)  errorBoundMode =  params->errorBoundMode;

        if(params->absErrBound >= 0) absErrBound = params->absErrBound;
        
        if(params->relBoundRatio >= 0) relBoundRatio = params->relBoundRatio;
        
        if(params->psnr >= 0) psnr = params->psnr;
        
		if(params->quantization_intervals>0)
		{
			updateQuantizationInfo(params->quantization_intervals);
			optQuantMode = 0;
		}
		else
			optQuantMode = 1;
	
		if(params->layers >= 0)
			layers = params->layers;
		if(params->sampleDistance >= 0)
			sampleDistance = params->sampleDistance;
		if(params->predThreshold > 0)
			predThreshold = params->predThreshold;
		if(params->psnr > 0)
			psnr = params->psnr;
		if(params->pw_relBoundRatio > 0)
			pw_relBoundRatio = params->pw_relBoundRatio;
		if(params->segment_size > 0)
			segment_size = params->segment_size;
		if(params->pwr_type >= 0)
			pwr_type = params->pwr_type;
    }

//	versionNumber[0] = SZ_VER_MAJOR; //0
//	versionNumber[1] = SZ_VER_MINOR; //5
//	versionNumber[2] = SZ_VER_REVISION; //15

	if(params->quantization_intervals%2!=0)
	{
		printf("Error: quantization_intervals must be an even number!\n");
		return SZ_NSCS;
	}

    //initialization for Huffman encoding
	SZ_Reset();
	
    return SZ_SCES;
}

int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	int dimension;
	if(r1==0) 
	{
		dimension = 0;
	}
	else if(r2==0) 
	{
		dimension = 1;
	}
	else if(r3==0) 
	{
		dimension = 2;
	}
	else if(r4==0) 
	{
		dimension = 3;
	}
	else if(r5==0) 
	{
		dimension = 4;
	}
	else 
	{
		dimension = 5;
	}
	return dimension;	
}

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t dataLength;
	if(r1==0) 
	{
		dataLength = 0;
	}
	else if(r2==0) 
	{
		dataLength = r1;
	}
	else if(r3==0) 
	{
		dataLength = r1*r2;
	}
	else if(r4==0) 
	{
		dataLength = r1*r2*r3;
	}
	else if(r5==0) 
	{
		dataLength = r1*r2*r3*r4;
	}
	else 
	{
		dataLength = r1*r2*r3*r4*r5;
	}
	return dataLength;
}

/*-------------------------------------------------------------------------*/
/**
    @brief      Perform Compression 
    @param      data           data to be compressed
    @param      outSize        the size (in bytes) after compression
    @param		r5,r4,r3,r2,r1	the sizes of each dimension (supporting only 5 dimensions at most in this version.
    @return     compressed data (in binary stream) or NULL(0) if any errors

 **/
/*-------------------------------------------------------------------------*/
unsigned char* SZ_compress_args(int dataType, void *data, size_t *outSize, int errBoundMode, double absErrBound, 
double relBoundRatio, double pwrBoundRatio, int pwrType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	//TODO
	conf_params->dataType = dataType;
	if(dataType==SZ_FLOAT)
	{
		unsigned char *newByteData = NULL;
		
		SZ_compress_args_float(&newByteData, (float *)data, r5, r4, r3, r2, r1, 
		outSize, errBoundMode, absErrBound, relBoundRatio, pwrBoundRatio, pwrType);
		
		return newByteData;
	}
	else if(dataType==SZ_DOUBLE)
	{
		unsigned char *newByteData;
		SZ_compress_args_double(&newByteData, (double *)data, r5, r4, r3, r2, r1, 
		outSize, errBoundMode, absErrBound, relBoundRatio, pwrBoundRatio, pwrType);
		
		return newByteData;
	}
	else if(dataType==SZ_INT64)
	{
		unsigned char *newByteData;
		SZ_compress_args_int64(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;
	}		
	else if(dataType==SZ_INT32) //int type
	{
		unsigned char *newByteData;
		SZ_compress_args_int32(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;
	}
	else if(dataType==SZ_INT16)
	{
		unsigned char *newByteData;
		SZ_compress_args_int16(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;		
	}
	else if(dataType==SZ_INT8)
	{
		unsigned char *newByteData;
		SZ_compress_args_int8(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;
	}
	else if(dataType==SZ_UINT64)
	{
		unsigned char *newByteData;
		SZ_compress_args_uint64(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;
	}		
	else if(dataType==SZ_UINT32) //int type
	{
		unsigned char *newByteData;
		SZ_compress_args_uint32(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;
	}
	else if(dataType==SZ_UINT16)
	{
		unsigned char *newByteData;
		SZ_compress_args_uint16(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;		
	}
	else if(dataType==SZ_UINT8)
	{
		unsigned char *newByteData;
		SZ_compress_args_uint8(&newByteData, data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
		return newByteData;
	} 	
	else
	{
		printf("Error: dataType can only be SZ_FLOAT, SZ_DOUBLE, SZ_INT8/16/32/64 or SZ_UINT8/16/32/64.\n");
		return NULL;
	}
}

int SZ_compress_args2(int dataType, void *data, unsigned char* compressed_bytes, size_t *outSize, 
int errBoundMode, double absErrBound, double relBoundRatio, double pwrBoundRatio, int pwrType, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	unsigned char* bytes = SZ_compress_args(dataType, data, outSize, errBoundMode, absErrBound, relBoundRatio, pwrBoundRatio, pwrType, r5, r4, r3, r2, r1);
    memcpy(compressed_bytes, bytes, *outSize);
    free(bytes); 
	return SZ_SCES;
}

int SZ_compress_args3(int dataType, void *data, unsigned char* compressed_bytes, size_t *outSize, int errBoundMode, double absErrBound, double relBoundRatio, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
size_t s5, size_t s4, size_t s3, size_t s2, size_t s1,
size_t e5, size_t e4, size_t e3, size_t e2, size_t e1)
{
	conf_params->dataType = dataType;
	if(dataType==SZ_FLOAT)
	{
		SZ_compress_args_float_subblock(compressed_bytes, (float *)data, 
		r5, r4, r3, r2, r1,
		s5, s4, s3, s2, s1,
		e5, e4, e3, e2, e1,
		outSize, errBoundMode, absErrBound, relBoundRatio);
		
		return SZ_SCES;
	}
	else if(dataType==SZ_DOUBLE)
	{
		SZ_compress_args_double_subblock(compressed_bytes, (double *)data, 
		r5, r4, r3, r2, r1,
		s5, s4, s3, s2, s1,
		e5, e4, e3, e2, e1,
		outSize, errBoundMode, absErrBound, relBoundRatio);
		
		return SZ_SCES;
	}
	else
	{
		printf("Error (in SZ_compress_args3): dataType can only be SZ_FLOAT or SZ_DOUBLE.\n");
		return SZ_NSCS;
	}	
}

unsigned char *SZ_compress(int dataType, void *data, size_t *outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{	
	unsigned char *newByteData = SZ_compress_args(dataType, data, outSize, errorBoundMode, absErrBound, relBoundRatio, 
	pw_relBoundRatio, pwr_type, r5, r4, r3, r2, r1);
	return newByteData;
}

//////////////////
/*-------------------------------------------------------------------------*/
/**
    @brief      Perform Compression 
    @param      data           data to be compressed
    @param		reservedValue  the reserved value
    @param      outSize        the size (in bytes) after compression
    @param		r5,r4,r3,r2,r1	the sizes of each dimension (supporting only 5 dimensions at most in this version.
    @return     compressed data (in binary stream)

 **/
/*-------------------------------------------------------------------------*/
unsigned char *SZ_compress_rev_args(int dataType, void *data, void *reservedValue, size_t *outSize, int errBoundMode, double absErrBound, double relBoundRatio, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t dataLength;
	unsigned char *newByteData;
	dataLength = computeDataLength(r5,r4,r3,r2,r1);
	//TODO
	printf("SZ compression with reserved data is TO BE DONE LATER.\n");
	exit(0);
	
	return newByteData;	
}

int SZ_compress_rev_args2(int dataType, void *data, void *reservedValue, unsigned char* compressed_bytes, size_t *outSize, int errBoundMode, double absErrBound, double relBoundRatio, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	conf_params->dataType = dataType;
	unsigned char* bytes = SZ_compress_rev_args(dataType, data, reservedValue, outSize, errBoundMode, absErrBound, relBoundRatio, r5, r4, r3, r2, r1);
	memcpy(compressed_bytes, bytes, *outSize);
	free(bytes); //free(bytes) is removed , because of dump error at MIRA system (PPC architecture), fixed?
	return 0;
}

unsigned char *SZ_compress_rev(int dataType, void *data, void *reservedValue, size_t *outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	size_t dataLength;

	unsigned char *newByteData;
	dataLength = computeDataLength(r5,r4,r3,r2,r1);	
	//TODO
	printf("SZ compression with reserved data is TO BE DONE LATER.\n");
	exit(0);
	
	return newByteData;
}

void *SZ_decompress(int dataType, unsigned char *bytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	int x = 1;
	char *y = (char*)&x;
	if(*y==1)
		sysEndianType = LITTLE_ENDIAN_SYSTEM;
	else //=0
		sysEndianType = BIG_ENDIAN_SYSTEM;
	
	if(dataType == SZ_FLOAT)
	{
		float *newFloatData;
		SZ_decompress_args_float(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
		return newFloatData;	
	}
	else if(dataType == SZ_DOUBLE)
	{
		double *newDoubleData;
		SZ_decompress_args_double(&newDoubleData, r5, r4, r3, r2, r1, bytes, byteLength);
		return newDoubleData;	
	}
	else if(dataType == SZ_INT8)
	{
		int8_t *newInt8Data;
		SZ_decompress_args_int8(&newInt8Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newInt8Data;
	}
	else if(dataType == SZ_INT16)
	{
		int16_t *newInt16Data;
		SZ_decompress_args_int16(&newInt16Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newInt16Data;
	}
	else if(dataType == SZ_INT32)
	{
		int32_t *newInt32Data;
		SZ_decompress_args_int32(&newInt32Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newInt32Data;
	}
	else if(dataType == SZ_INT64)
	{
		int64_t *newInt64Data;
		SZ_decompress_args_int64(&newInt64Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newInt64Data;
	}
	else if(dataType == SZ_UINT8)
	{
		uint8_t *newUInt8Data;
		SZ_decompress_args_uint8(&newUInt8Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newUInt8Data;
	}
	else if(dataType == SZ_UINT16)
	{
		uint16_t *newUInt16Data;
		SZ_decompress_args_uint16(&newUInt16Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newUInt16Data;
	}
	else if(dataType == SZ_UINT32)
	{
		uint32_t *newUInt32Data;
		SZ_decompress_args_uint32(&newUInt32Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newUInt32Data;
	}
	else if(dataType == SZ_UINT64)
	{
		uint64_t *newUInt64Data;
		SZ_decompress_args_uint64(&newUInt64Data, r5, r4, r3, r2, r1, bytes, byteLength);
		return newUInt64Data;
	}
	else 
	{
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		return NULL;	
	}
}

/**
 * 
 * 
 * return number of elements or -1 if any errors
 * */
size_t SZ_decompress_args(int dataType, unsigned char *bytes, size_t byteLength, void* decompressed_array, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	//size_t i;
	size_t nbEle = computeDataLength(r5,r4,r3,r2,r1);
	
	if(dataType == SZ_FLOAT)
	{
		float* data = (float *)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		float* data_array = (float *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(float));
		//for(i=0;i<nbEle;i++)
		//	data_array[i] = data[i];	
		free(data); //this free operation seems to not work with BlueG/Q system.	
	}
	else if (dataType == SZ_DOUBLE)
	{
		double* data = (double *)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		double* data_array = (double *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(double));
		//for(i=0;i<nbEle;i++)
		//	data_array[i] = data[i];
		free(data); //this free operation seems to not work with BlueG/Q system.	
	}
	else if(dataType == SZ_INT8)
	{
		int8_t* data = (int8_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		int8_t* data_array = (int8_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(int8_t));
		free(data);
	}
	else if(dataType == SZ_INT16)
	{
		int16_t* data = (int16_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		int16_t* data_array = (int16_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(int16_t));
		free(data);	
	}
	else if(dataType == SZ_INT32)
	{
		int32_t* data = (int32_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		int32_t* data_array = (int32_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(int32_t));
		free(data);	
	}
	else if(dataType == SZ_INT64)
	{
		int64_t* data = (int64_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		int64_t* data_array = (int64_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(int64_t));
		free(data);		
	}
	else if(dataType == SZ_UINT8)
	{
		uint8_t* data = (uint8_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		uint8_t* data_array = (uint8_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(uint8_t));
		free(data);
	}
	else if(dataType == SZ_UINT16)
	{
		uint16_t* data = (uint16_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		uint16_t* data_array = (uint16_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(uint16_t));
		free(data);		
	}
	else if(dataType == SZ_UINT32)
	{
		uint32_t* data = (uint32_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		uint32_t* data_array = (uint32_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(uint32_t));
		free(data);		
	}
	else if(dataType == SZ_UINT64)
	{
		uint64_t* data = (uint64_t*)SZ_decompress(dataType, bytes, byteLength, r5, r4, r3, r2, r1);
		uint64_t* data_array = (uint64_t *)decompressed_array;
		memcpy(data_array, data, nbEle*sizeof(uint64_t));
		free(data);			
	}
	else
	{ 
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		return SZ_NSCS; //indicating error		
	}

	return nbEle;
}

sz_metadata* SZ_getMetadata(unsigned char* bytes)
{
	int index = 0, i, isConstant, isLossless, sizeType;
	size_t dataSeriesLength = 0;
	int versions[3] = {0,0,0};
	for (i = 0; i < 3; i++)
		versions[i] = bytes[index++]; //3
	unsigned char sameRByte = bytes[index++]; //1
	isConstant = sameRByte & 0x01;
	//szMode = (sameRByte & 0x06)>>1;
	isLossless = (sameRByte & 0x10)>>4;
	SZ_SIZE_TYPE = sizeType = ((sameRByte & 0x40)>>6)==1?8:4;
	
	sz_params* params = convertBytesToSZParams(&(bytes[index]));
	if(conf_params!=NULL)
		free(conf_params);
	conf_params = params;	
	index += MetaDataByteLength;
	
	if(params->dataType!=SZ_FLOAT && params->dataType!= SZ_DOUBLE) //if this type is an Int type
		index++; //jump to the dataLength info byte address
	dataSeriesLength = bytesToSize(&(bytes[index]));// 4 or 8	
	
	sz_metadata* metadata = (sz_metadata*)malloc(sizeof(struct sz_metadata));
	
	metadata->versionNumber[0] = versions[0];
	metadata->versionNumber[1] = versions[1];
	metadata->versionNumber[2] = versions[2];
	metadata->isConstant = isConstant;
	metadata->isLossless = isLossless;
	metadata->sizeType = sizeType;
	metadata->dataSeriesLength = dataSeriesLength;
	
	metadata->conf_params = conf_params;
	
	return metadata;
}

void SZ_printMetadata(sz_metadata* metadata)
{
	printf("=================SZ Compression Meta Data=================\n");
	printf("Version:                        \t %d.%d.%d\n", metadata->versionNumber[0], metadata->versionNumber[1], metadata->versionNumber[2]);
	printf("Constant data?:                 \t %s\n", metadata->isConstant==1?"YES":"NO");
	printf("Lossless?:                      \t %s\n", metadata->isLossless==1?"YES":"NO");
	printf("Size type (size of # elements): \t %d bytes\n", metadata->sizeType); 
	printf("Num of elements:                \t %zu\n", metadata->dataSeriesLength);
		
	sz_params* params = metadata->conf_params;
	
	switch(params->dataType)
	{
	case SZ_FLOAT:
		printf("Data type:                      \t FLOAT\n");
		break;
	case SZ_DOUBLE:
		printf("Data type:                      \t DOUBLE\n");
		break;
	case SZ_INT8:
		printf("Data type:                      \t INT8\n");
		break;	
	case SZ_INT16:
		printf("Data type:                      \t INT16\n");
		break;
	case SZ_INT32:
		printf("Data type:                      \t INT32\n");
		break;	
	case SZ_INT64:
		printf("Data type:                      \t INT64\n");
		break;	
	case SZ_UINT8:
		printf("Data type:                      \t UINT8\n");
		break;	
	case SZ_UINT16:
		printf("Data type:                      \t UINT16\n");
		break;
	case SZ_UINT32:
		printf("Data type:                      \t UINT32\n");
		break;	
	case SZ_UINT64:
		printf("Data type:                      \t UINT64\n");
		break;				
	}
	
	if(optQuantMode==1)
	{
		printf("quantization_intervals:         \t 0\n");
		printf("max_quant_intervals:            \t %d\n", params->max_quant_intervals);
	}
	else
	{
		printf("quantization_intervals:         \t %d\n", params->quantization_intervals);
		printf("max_quant_intervals:            \t - %d\n");		
	}
	
	printf("dataEndianType (prior raw data):\t %s\n", params->dataEndianType==1?"BIG_ENDIAN":"LITTLE_ENDIAN");
	printf("sysEndianType (at compression): \t %s\n", params->sysEndianType==1?"BIG_ENDIAN":"LITTLE_ENDIAN");
	printf("sampleDistance:                 \t %d\n", params->sampleDistance);
	printf("predThreshold:                  \t %f\n", params->predThreshold);
	switch(params->szMode)
	{
	case SZ_BEST_SPEED:
		printf("szMode:                         \t SZ_BEST_SPEED (without Gzip)\n");
		break;
	case SZ_BEST_COMPRESSION:
		printf("szMode:                         \t SZ_BEST_COMPRESSION (with Gzip)\n");
		break;
	}
	switch(params->gzipMode)
	{
	case Z_BEST_SPEED:
		printf("gzipMode:                       \t Z_BEST_SPEED\n");
		break;
	case Z_DEFAULT_COMPRESSION:
		printf("gzipMode:                       \t Z_BEST_SPEED\n");
		break;	
	case Z_BEST_COMPRESSION:
		printf("gzipMode:                       \t Z_BEST_COMPRESSION\n");
		break;
	}
	
	switch(params->errorBoundMode)
	{
	case ABS:
		printf("errBoundMode:                   \t ABS\n");
		printf("absErrBound:                    \t %f\n", params->absErrBound);
		break;
	case REL:
		printf("errBoundMode:                   \t REL (based on value_range extent)\n");
		printf("relBoundRatio:                  \t %f\n", params->relBoundRatio);
		break;
	case ABS_AND_REL:
		printf("errBoundMode:                   \t ABS_AND_REL\n");
		printf("absErrBound:                    \t %f\n", params->absErrBound);
		printf("relBoundRatio:                  \t %f\n", params->relBoundRatio);
		break;
	case ABS_OR_REL:
		printf("errBoundMode:                   \t ABS_OR_REL\n");
		printf("absErrBound:                    \t %f\n", params->absErrBound);
		printf("relBoundRatio:                  \t %f\n", params->relBoundRatio);
		break;
	case PSNR:
		printf("errBoundMode:                   \t PSNR\n");
		printf("psnr:                           \t %f\n", params->psnr);
		break;
	case PW_REL:
		printf("errBoundMode:                   \t PW_REL\n");
		break;
	case ABS_AND_PW_REL:
		printf("errBoundMode:                   \t ABS_AND_PW_REL\n");
		printf("absErrBound:                    \t %f\n", params->absErrBound);
		break;
	case ABS_OR_PW_REL:
		printf("errBoundMode:                   \t ABS_OR_PW_REL\n");
		printf("absErrBound:                    \t %f\n", params->absErrBound);
		break;
	case REL_AND_PW_REL:
		printf("errBoundMode:                   \t REL_AND_PW_REL\n");
		printf("range_relBoundRatio:            \t %f\n", params->relBoundRatio);
		break;
	case REL_OR_PW_REL:
		printf("errBoundMode:                   \t REL_OR_PW_REL\n");
		printf("range_relBoundRatio:            \t %f\n", params->relBoundRatio);
		break;
	}
	
	if(params->errorBoundMode>=PW_REL && params->errorBoundMode<=REL_OR_PW_REL)
	{
		printf("pw_relBoundRatio:               \t %f\n", params->pw_relBoundRatio);
		printf("segment_size:                   \t %d\n", params->segment_size);
		switch(params->pwr_type)
		{
		case SZ_PWR_MIN_TYPE:
			printf("pwrType:                    \t SZ_PWR_MIN_TYPE\n");
			break;
		case SZ_PWR_AVG_TYPE:
			printf("pwrType:                    \t SZ_PWR_AVG_TYPE\n");
			break;
		case SZ_PWR_MAX_TYPE:
			printf("pwrType:                    \t SZ_PWR_MAX_TYPE\n");
			break;
		}
	}
}

/*-----------------------------------batch data compression--------------------------------------*/

void filloutDimArray(size_t* dim, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	if(r2==0)
		dim[0] = r1;
	else if(r3==0)
	{
		dim[0] = r2;
		dim[1] = r1;
	}
	else if(r4==0)
	{
		dim[0] = r3;
		dim[1] = r2;
		dim[2] = r1;
	}
	else if(r5==0)
	{
		dim[0] = r4;
		dim[1] = r3;
		dim[2] = r2;
		dim[3] = r1;
	}
	else
	{
		dim[0] = r5;
		dim[1] = r4;
		dim[2] = r3;
		dim[3] = r2;
		dim[4] = r1;		
	}
}

size_t compute_total_batch_size()
{
	size_t eleNum = 0, totalSize = 0;
	SZ_Variable* p = sz_varset->header;
	while(p->next!=NULL)
	{
		eleNum = computeDataLength(p->next->r5, p->next->r4, p->next->r3, p->next->r2, p->next->r1);
		if(p->next->dataType==SZ_FLOAT)
			totalSize += (eleNum*4);
		else
			totalSize += (eleNum*8);
		p=p->next;
	}
	return totalSize;
}

int isZlibFormat(unsigned char magic1, unsigned char magic2)
{
	if(magic1==104&&magic2==5) //DC+BS
		return 1;
	if(magic1==104&&magic2==129) //DC+DC
		return 1;
	if(magic1==104&&magic2==222) //DC+BC
		return 1;
	if(magic1==120&&magic2==1) //BC+BS
		return 1;
	if(magic1==120&&magic2==156) //BC+DC
		return 1;
	if(magic1==120&&magic2==218) //BC+BS
		return 1;
	return 0;
}

unsigned char* SZ_batch_compress(size_t *outSize)
{	
	size_t dataLength;
	DynamicByteArray* dba; 
	new_DBA(&dba, 32768);
	
	//number of variables
	int varCount = sz_varset->count;
	unsigned char countBufBytes[SZ_SIZE_TYPE];
	intToBytes_bigEndian(countBufBytes, varCount);
	//add only the lats two bytes, i.e., the maximum # variables is supposed to be less than 32768
	addDBA_Data(dba, countBufBytes[2]);
	addDBA_Data(dba, countBufBytes[3]);
	
	size_t i, j, k = 0;
	SZ_Variable* p = sz_varset->header->next;
	while(p!=NULL)
	{
		if(p->dataType==SZ_FLOAT)
		{
			unsigned char *newByteData;
			size_t outSize;
			SZ_compress_args_float_wRngeNoGzip(&newByteData, (float *)p->data, 
			p->r5, p->r4, p->r3, p->r2, p->r1, 
			&(p->compressedSize), p->errBoundMode, p->absErrBound, p->relBoundRatio, p->pwRelBoundRatio);
			
			p->compressedBytes = newByteData;
		}
		else if(p->dataType==SZ_DOUBLE)
		{
			unsigned char *newByteData;
			size_t outSize;
			SZ_compress_args_double_wRngeNoGzip(&newByteData, (double *)p->data, 
			p->r5, p->r4, p->r3, p->r2, p->r1, 
			&(p->compressedSize), p->errBoundMode, p->absErrBound, p->relBoundRatio, p->pwRelBoundRatio);
			
			p->compressedBytes = newByteData;
		}
		
		SZ_ReleaseHuffman();
		
		//TODO: copy metadata information (totally 1+4*x+4+(1+y) bytes)
		//byte format of variable storage: meta 1 bytes (000000xy) x=0 SZ/1 HZ, y=0 float/1 double ; 
		//4*x: x integers indicating the dimension sizes; 
		//compressed length (int) 4 bytes;
		//1+y: 1 means the length of varname, y bytes represents the varName string. 
		int meta = 0;
		if(p->dataType == SZ_FLOAT)
			meta = 2; //10
		else if(p->dataType == SZ_DOUBLE)
			meta = 3; //11
		
		//keep dimension information
		int dimNum = computeDimension(p->r5, p->r4, p->r3, p->r2, p->r1);
		size_t dimSize[dimNum];
		memset(dimSize, 0, dimNum*SZ_SIZE_TYPE);
		meta = meta | dimNum << 2; //---aaabc: aaa indicates dim, b indicates HZ, c indicates dataType
		
		addDBA_Data(dba, (unsigned char)meta);
		
		filloutDimArray(dimSize, p->r5, p->r4, p->r3, p->r2, p->r1);
		
		for(j=0;j<dimNum;j++)
		{
			if(SZ_SIZE_TYPE==4)
				intToBytes_bigEndian(countBufBytes, dimSize[j]);
			else
				longToBytes_bigEndian(countBufBytes, dimSize[j]);
			for(i = 0;i<SZ_SIZE_TYPE;i++)
				addDBA_Data(dba, countBufBytes[i]);
		}
			 
		//Keep compressed size information	 
		if(SZ_SIZE_TYPE==4)
			intToBytes_bigEndian(countBufBytes, p->compressedSize);
		else
			longToBytes_bigEndian(countBufBytes, p->compressedSize);
		
		for(i = 0;i<SZ_SIZE_TYPE;i++)
			addDBA_Data(dba, countBufBytes[i]);			 
			 
		//Keep varName information	 
		unsigned char varNameLength = (unsigned char)strlen(p->varName);
		addDBA_Data(dba, varNameLength);
		memcpyDBA_Data(dba, (unsigned char*)p->varName, varNameLength);
			 
		//copy the compressed stream
		memcpyDBA_Data(dba, p->compressedBytes, p->compressedSize);
		free(p->compressedBytes);
		p->compressedBytes = NULL;
			 
		p = p->next;
	}
	
	unsigned char* tmpFinalCompressedBytes;
	convertDBAtoBytes(dba, &tmpFinalCompressedBytes);
	
	unsigned char* tmpCompressedBytes2;
	size_t tmpGzipSize = 0;
	
	if(szMode!=SZ_BEST_SPEED)
	{
		tmpGzipSize = zlib_compress5(tmpFinalCompressedBytes, dba->size, &tmpCompressedBytes2, gzipMode);
		free(tmpFinalCompressedBytes);		
	}
	else
	{
		tmpCompressedBytes2 = tmpFinalCompressedBytes;
		tmpGzipSize = dba->size;
	}	
	
	unsigned char* finalCompressedBytes = (unsigned char*) malloc(sizeof(unsigned char)*(SZ_SIZE_TYPE+tmpGzipSize));
	
	if(SZ_SIZE_TYPE==4)
		intToBytes_bigEndian(countBufBytes, dba->size);
	else
		longToBytes_bigEndian(countBufBytes, dba->size);
	
	memcpy(finalCompressedBytes, countBufBytes, SZ_SIZE_TYPE);
	
	memcpy(&(finalCompressedBytes[SZ_SIZE_TYPE]), tmpCompressedBytes2, tmpGzipSize);
	free(tmpCompressedBytes2);
	
	*outSize = SZ_SIZE_TYPE+tmpGzipSize;
	free_DBA(dba);
	
	return finalCompressedBytes;
}

SZ_VarSet* SZ_batch_decompress(unsigned char* compressedStream, size_t compressedLength, int *status)
{
	size_t i, j, k = 0;
	unsigned char sizeByteBuf[SZ_SIZE_TYPE];
	
	int x = 1;
	char *y = (char*)&x;
	if(*y==1)
		sysEndianType = LITTLE_ENDIAN_SYSTEM;
	else //=0
		sysEndianType = BIG_ENDIAN_SYSTEM;
	
	//get target decompression size for Gzip (zlib)
	sizeByteBuf[0] = compressedStream[0];
	sizeByteBuf[1] = compressedStream[1];
	sizeByteBuf[2] = compressedStream[2];
	sizeByteBuf[3] = compressedStream[3];
	if(SZ_SIZE_TYPE==8)
	{
		sizeByteBuf[4] = compressedStream[4];
		sizeByteBuf[5] = compressedStream[5];
		sizeByteBuf[6] = compressedStream[6];
		sizeByteBuf[7] = compressedStream[7];		
	}
	
	size_t targetUncompressSize = 0; 
	if(SZ_SIZE_TYPE==4)
		targetUncompressSize = bytesToInt_bigEndian(sizeByteBuf);
	else
		targetUncompressSize = bytesToLong_bigEndian(sizeByteBuf);
	
	//Gzip decompression
	unsigned char* gzipDecpressBytes;
	size_t gzipDecpressSize = 0;
	if(isZlibFormat(compressedStream[SZ_SIZE_TYPE], compressedStream[SZ_SIZE_TYPE+1])!=0)
	{
		gzipDecpressSize = zlib_uncompress5(&(compressedStream[SZ_SIZE_TYPE]), (unsigned long)compressedLength, &gzipDecpressBytes, (unsigned long)targetUncompressSize);
	}
	else
	{
		gzipDecpressSize = compressedLength;
		gzipDecpressBytes = &(compressedStream[SZ_SIZE_TYPE]);
	}
	
	if(gzipDecpressSize!=targetUncompressSize)
	{
		printf("Error: incorrect decompression in zlib_uncompress3: gzipDecpressSize!=targetUncompressSize\n");
		*status = SZ_NSCS;
		return NULL;
	}
	
	//Start analyzing the byte stream for further decompression	
	sizeByteBuf[0] = 0;
	sizeByteBuf[1] = 0; 
	sizeByteBuf[2] = gzipDecpressBytes[k++];
	sizeByteBuf[3] = gzipDecpressBytes[k++];
	
	int varCount = bytesToInt_bigEndian(sizeByteBuf);	
		
	int varNum = sz_varset->count;
	size_t dataLength, cpressedLength;
	
	if(varNum==0)
	{
		SZ_Variable* lastVar = sz_varset->header; 
		if(lastVar->next!=NULL)
		{
			printf("Error: sz_varset.count is inconsistent with the number of variables in sz_varset->header.\n");
			*status = SZ_NSCS;
			return NULL;
		}
		for(i=0;i<varCount;i++)
		{
			int type = (int)gzipDecpressBytes[k++];
			int dataType;
			int tmpType = type & 0x03;
			if(tmpType==2) //FLOAT
				dataType = SZ_FLOAT;
			else if(tmpType==3)//DOUBLE
				dataType = SZ_DOUBLE;
			else
			{
				printf("Error: Wrong value of decompressed type. \n Please check the correctness of the decompressed data.\n");
				*status = SZ_NSCS;
				return NULL;
			}
			
			//get # dimensions and the size of each dimension
			int dimNum = (type & 0x1C) >> 2; //compute dimension
			int start_dim = 5 - dimNum;
			size_t dimSize[5];
			memset(dimSize, 0, 5*SZ_SIZE_TYPE);
			
			for(j=0;j<dimNum;j++)
			{
				memcpy(sizeByteBuf, &(gzipDecpressBytes[k]), SZ_SIZE_TYPE);
				k+=SZ_SIZE_TYPE;
				if(SZ_SIZE_TYPE==4)
					dimSize[start_dim+j] = bytesToInt_bigEndian(sizeByteBuf);
				else
					dimSize[start_dim+j] = bytesToLong_bigEndian(sizeByteBuf);
			}	
			
			//get compressed length
			memcpy(sizeByteBuf, &(gzipDecpressBytes[k]), SZ_SIZE_TYPE);
			k+=SZ_SIZE_TYPE;
			if(SZ_SIZE_TYPE==4)
				cpressedLength = bytesToInt_bigEndian(sizeByteBuf);	
			else
				cpressedLength = bytesToLong_bigEndian(sizeByteBuf);	
					
			//Keep varName information	 
			int varNameLength = gzipDecpressBytes[k++];
			char* varNameString = (char*)malloc(sizeof(char)*(varNameLength+1));	
			memcpy(varNameString, &(gzipDecpressBytes[k]), varNameLength);
			k+=varNameLength;
			varNameString[varNameLength] = '\0';
		
			//TODO: convert szTmpBytes to data array.
			dataLength = computeDataLength(dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4]);
			int dim = computeDimension(dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4]);
			if(dataType==SZ_FLOAT)
			{
				float* newData;
				TightDataPointStorageF* tdps;
				int errBoundMode = new_TightDataPointStorageF_fromFlatBytes(&tdps, &(gzipDecpressBytes[k]), cpressedLength);

				if (dim == 1)
					getSnapshotData_float_1D(&newData,dimSize[4],tdps, errBoundMode);
				else
				if (dim == 2)
					getSnapshotData_float_2D(&newData,dimSize[3],dimSize[4],tdps, errBoundMode);
				else
				if (dim == 3)
					getSnapshotData_float_3D(&newData,dimSize[2], dimSize[3], dimSize[4],tdps, errBoundMode);
				else
				if (dim == 4)
					getSnapshotData_float_3D(&newData,dimSize[1]*dimSize[2], dimSize[3], dimSize[4],tdps, errBoundMode);
				else
				{
					printf("Current version doesn't support 5 dimensions.\n");
					*status = SZ_DERR;
					return NULL;
				}
				
				SZ_batchAddVar(varNameString, dataType, newData, 0, 0, 0, dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4]);
				free_TightDataPointStorageF(tdps);			
			}	
			else if(dataType==SZ_DOUBLE)
			{
				double* newData;
				TightDataPointStorageD* tdps;
				int errBoundMode = new_TightDataPointStorageD_fromFlatBytes(&tdps, &(gzipDecpressBytes[k]), cpressedLength);

				if (dim == 1)
					getSnapshotData_double_1D(&newData,dimSize[4],tdps, errBoundMode);
				else
				if (dim == 2)
					getSnapshotData_double_2D(&newData,dimSize[3],dimSize[4],tdps, errBoundMode);
				else
				if (dim == 3)
					getSnapshotData_double_3D(&newData,dimSize[2], dimSize[3], dimSize[4],tdps, errBoundMode);
				else
				if (dim == 4)
					getSnapshotData_double_3D(&newData,dimSize[1]*dimSize[2], dimSize[3], dimSize[4],tdps, errBoundMode);
				else
				{
					printf("Current version doesn't support 5 dimensions.\n");
					*status = SZ_DERR;
					return NULL;
				}
			
				SZ_batchAddVar(varNameString, dataType, newData, 0, 0, 0, dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4]);
				free_TightDataPointStorageD(tdps);					
			}
			else
			{
				printf("Error: wrong dataType in the batch decompression.\n");
				*status = SZ_TERR;
				return NULL;
			}
	
			k+=cpressedLength;
		}
	}
	else if(varNum!=varCount)
	{
		printf("Error: Inconsistency of the # variables between sz_varset and decompressed stream.\n");
		printf("Note sz_varset.count should be either 0 or the correct number of variables stored in the decompressed stream.\n");
		printf("Currently, sz_varset.count=%d, expected number of variables = %d\n", varNum, varCount);
		*status = SZ_NSCS;
		return NULL;
	}
	else //if(varNum>0 && varNum==varCount)
	{
		SZ_Variable* p = sz_varset->header; 
		for(i=0;i<varCount;i++)
		{
			int type = (int)gzipDecpressBytes[k++];
			int dataType;
			int tmpType = type & 0x03;
			if(tmpType==2) //FLOAT
				dataType = SZ_FLOAT;
			else if(tmpType==3)//DOUBLE
				dataType = SZ_DOUBLE;
			else
			{
				printf("Error: Wrong value of decompressed type. \n Please check the correctness of the decompressed data.\n");
				*status = SZ_DERR;
				return NULL;
			}
			
			//get # dimensions and the size of each dimension
			int dimNum = (type & 0x1C) >> 2; //compute dimension
			int start_dim = 5 - dimNum;
			size_t dimSize[5];
			memset(dimSize, 0, 5*SZ_SIZE_TYPE);
			
			for(j=0;j<dimNum;j++)
			{
				memcpy(sizeByteBuf, &(gzipDecpressBytes[k]), SZ_SIZE_TYPE);
				k += SZ_SIZE_TYPE;
				if(SZ_SIZE_TYPE==4)
					dimSize[start_dim+j] = bytesToInt_bigEndian(sizeByteBuf);
				else
					dimSize[start_dim+j] = bytesToLong_bigEndian(sizeByteBuf);
			}
			
			//get compressed length
			memcpy(sizeByteBuf, &(gzipDecpressBytes[k]), SZ_SIZE_TYPE);
			k+=SZ_SIZE_TYPE;
			if(SZ_SIZE_TYPE==4)
				cpressedLength = bytesToInt_bigEndian(sizeByteBuf);	
			else
				cpressedLength = bytesToLong_bigEndian(sizeByteBuf);	
			
			//Keep varName information	 
			int varNameLength = gzipDecpressBytes[k++];
			char* varNameString = (char*)malloc(sizeof(char)*(varNameLength+1));	
			memcpy(varNameString, &(gzipDecpressBytes[k]), varNameLength);
			k+=varNameLength;
			varNameString[varNameLength] = '\0';
			
			SZ_Variable* curVar = p->next;
			
			int checkVarName = strcmp(varNameString, curVar->varName);
			if(checkVarName!=0)
			{
				printf("Error: the varNames in the compressed stream are inconsistent with those in the sz_varset\n");
				*status = SZ_DERR;
				return NULL;
			}
						
			//TODO: convert szTmpBytes to data array.
			dataLength = computeDataLength(dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4]);
			int dim = computeDimension(dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4]);			
			if(dataType==SZ_FLOAT)
			{
				float* newData;
				TightDataPointStorageF* tdps;
				int errBoundMode = new_TightDataPointStorageF_fromFlatBytes(&tdps, &(gzipDecpressBytes[k]), cpressedLength);					

				if (dim == 1)
					getSnapshotData_float_1D(&newData,dimSize[4],tdps, errBoundMode);
				else
				if (dim == 2)
					getSnapshotData_float_2D(&newData,dimSize[3],dimSize[4],tdps, errBoundMode);
				else
				if (dim == 3)
					getSnapshotData_float_3D(&newData,dimSize[2], dimSize[3], dimSize[4],tdps, errBoundMode);
				else
				if (dim == 4)
					getSnapshotData_float_3D(&newData,dimSize[1]*dimSize[2], dimSize[3], dimSize[4],tdps, errBoundMode);
				else
				{
					printf("Error: doesn't support 5 dimensions yet.\n");			
					*status = SZ_DERR;
					return NULL;
				}
				
				free_TightDataPointStorageF(tdps);	
				curVar->data = newData;
			}	
			else if(dataType==SZ_DOUBLE)
			{
				double* newData;
				TightDataPointStorageD* tdps;
				int errBoundMode = new_TightDataPointStorageD_fromFlatBytes(&tdps, &(gzipDecpressBytes[k]), cpressedLength);

				if (dim == 1)
					getSnapshotData_double_1D(&newData,dimSize[4],tdps, errBoundMode);
				else
				if (dim == 2)
					getSnapshotData_double_2D(&newData,dimSize[3],dimSize[4],tdps, errBoundMode);
				else
				if (dim == 3)
					getSnapshotData_double_3D(&newData,dimSize[2], dimSize[3], dimSize[4],tdps, errBoundMode);
				else
				if (dim == 4)
					getSnapshotData_double_3D(&newData,dimSize[1]*dimSize[2], dimSize[3], dimSize[4],tdps, errBoundMode);
				else
				{
					printf("Error: doesn't support 5 dimensions yet.\n");
					*status = SZ_DERR;
					return NULL;
				}
				SZ_batchAddVar(varNameString, dataType, newData, dimSize[0], dimSize[1], dimSize[2], dimSize[3], dimSize[4], 0, 0, 0);
				free_TightDataPointStorageD(tdps);	
			}
			else
			{
				printf("Error: wrong dataType in the batch decompression.\n");
				*status = SZ_TERR;
				return NULL;
			}		
			
			free(varNameString);
			k+=cpressedLength;
			p = p->next;
		}
	}
	free(gzipDecpressBytes);
	SZ_ReleaseHuffman();
	*status = SZ_SCES;
	return sz_varset;
}

void SZ_Finalize()
{
	//if(sz_varset!=NULL)
	//	free_VarSet();
	
	SZ_ReleaseHuffman();
	if(conf_params!=NULL)
		free(conf_params);
}

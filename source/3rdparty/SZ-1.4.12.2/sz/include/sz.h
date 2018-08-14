/**
 *  @file sz.h
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief Header file for the whole detector.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _SZ_H
#define _SZ_H

#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>      /* For gettimeofday(), in microseconds */
#include <time.h>          /* For time(), in seconds */
#include "iniparser.h"
#include "CompressElement.h"
#include "DynamicByteArray.h"
#include "DynamicIntArray.h"
#include "VarSet.h"
#include "Huffman.h"
#include "TightDataPointStorageD.h"
#include "TightDataPointStorageF.h"
#include "TightDataPointStorageI.h"
#include "conf.h"
#include "dataCompression.h"
#include "ByteToolkit.h"
#include "TypeManager.h"
#include "sz_int8.h"
#include "sz_int16.h"
#include "sz_int32.h"
#include "sz_int64.h"
#include "sz_uint8.h"
#include "sz_uint16.h"
#include "sz_uint32.h"
#include "sz_uint64.h"
#include "sz_float.h"
#include "sz_double.h"
#include "szd_int8.h"
#include "szd_int16.h"
#include "szd_int32.h"
#include "szd_int64.h"
#include "szd_uint8.h"
#include "szd_uint16.h"
#include "szd_uint32.h"
#include "szd_uint64.h"
#include "szd_float.h"
#include "szd_double.h"
#include "sz_float_pwr.h"
#include "sz_double_pwr.h"
#include "callZlib.h"
#include "rw.h"

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

//typedef char int8_t;
//typedef unsigned char uint8_t;
//typedef short int16_t;
//typedef unsigned short uint16_t;
//typedef int int32_t;
//typedef unsigned int uint32_t;
//typedef long int64_t;
//typedef unsigned long uint64_t;

#define SZ_VERNUM 0x0140
#define SZ_VER_MAJOR 1
#define SZ_VER_MINOR 4
#define SZ_VER_BUILD 12
#define SZ_VER_REVISION 1

#define HZ 102
#define SZ 101

#define ABS 0
#define REL 1
#define ABS_AND_REL 2
#define ABS_OR_REL 3
#define PSNR 4

#define PW_REL 10
#define ABS_AND_PW_REL 11
#define ABS_OR_PW_REL 12
#define REL_AND_PW_REL 13
#define REL_OR_PW_REL 14

#define SZ_FLOAT 0
#define SZ_DOUBLE 1
#define SZ_UINT8 2
#define SZ_INT8 3
#define SZ_UINT16 4
#define SZ_INT16 5
#define SZ_UINT32 6
#define SZ_INT32 7
#define SZ_UINT64 8
#define SZ_INT64 9

#define LITTLE_ENDIAN_DATA 0
#define BIG_ENDIAN_DATA 1 //big_endian (ppc, max, etc.) ; little_endian (x86, x64, etc.)

#define LITTLE_ENDIAN_SYSTEM 0
#define BIG_ENDIAN_SYSTEM 1

#define DynArrayInitLen 1024

#define MIN_ZLIB_DEC_ALLOMEM_BYTES 1000000

//#define maxRangeRadius 32768
//#define maxRangeRadius 1048576//131072

#define SZ_BEST_SPEED 0
#define SZ_BEST_COMPRESSION 1
#define SZ_DEFAULT_COMPRESSION 2

#define SZ_PWR_MIN_TYPE 0
#define SZ_PWR_AVG_TYPE 1
#define SZ_PWR_MAX_TYPE 2

//SUCCESS returning status
#define SZ_SCES 0  //successful
#define SZ_NSCS -1 //Not successful
#define SZ_FERR -2 //Failed to open input file
#define SZ_TERR -3 //wrong data type (should be only float or double)
#define SZ_DERR -4 //dimension error
#define SZ_MERR -5 //sz_mode error
#define SZ_BERR -6 //bound-mode error (should be only ABS, REL, ABS_AND_REL, ABS_OR_REL, or PW_REL)

#define SZ_MAINTAIN_VAR_DATA 0
#define SZ_DESTROY_WHOLE_VARSET 1

#define GROUP_COUNT 16 //2^{16}=65536
	
#define MetaDataByteLength 20	
	
//Note: the following setting should be consistent with stateNum in Huffman.h
//#define intvCapacity 65536
//#define intvRadius 32768
//#define intvCapacity 131072
//#define intvRadius 65536

extern unsigned int maxRangeRadius;

extern int intvCapacity;
extern int intvRadius;

extern int sysEndianType; //endian type of the system
extern int dataEndianType; //endian type of the data
//extern int maxSegmentNum;

extern char maxHeap[10];
 
extern long status;

extern int sol_ID;
extern int errorBoundMode; //ABS, REL, ABS_AND_REL, or ABS_OR_REL, PW_REL

extern int gzipMode; //four options: Z_NO_COMPRESSION, or Z_BEST_SPEED, Z_BEST_COMPRESSION, Z_DEFAULT_COMPRESSION

extern char *sz_cfgFile;

extern int offset;

extern double absErrBound;
extern double relBoundRatio;
extern double psnr;
extern double pw_relBoundRatio;
extern int segment_size;

extern int pwr_type;

extern int versionNumber[4];

extern int layers;
extern float predThreshold;
extern int sampleDistance;
extern char optQuantMode;

extern int szMode; //0 (best speed) or 1 (better compression with Gzip)

//extern int spaceFillingCurveTransform; //default is 0, or 1 set by sz.config
//extern int reOrgSize; //the granularity of the reganization of the original data

extern SZ_VarSet* sz_varset;

extern int SZ_SIZE_TYPE; //4 or 8: sizeof(size_t) 

//typedef unsigned long unsigned long;
//typedef unsigned int uint;

typedef union lint16
{
	unsigned short usvalue;
	short svalue;
	unsigned char byte[2];
} lint16;

typedef union lint32
{
	int ivalue;
	unsigned int uivalue;
	unsigned char byte[4];
} lint32;

typedef union lint64
{
	long lvalue;
	unsigned long ulvalue;
	unsigned char byte[8];
} lint64;

typedef union ldouble
{
    double value;
    unsigned long lvalue;
    unsigned char byte[8];
} ldouble;

typedef union lfloat
{
    float value;
    unsigned int ivalue;
    unsigned char byte[4];
} lfloat;

/* array meta data and compression parameters for SZ_Init_Params() */
typedef struct sz_params
{
	int dataType;
	unsigned int max_quant_intervals;
	unsigned int quantization_intervals;
    int dataEndianType; //*
    int sysEndianType; //*sysEndianType is actually set automatically.
    int sol_ID;//x
    int layers;//x
    int sampleDistance; //2 bytes
    float predThreshold;  // 2 bytes
    int offset;//x
    int szMode; //*
    int gzipMode; //*
    int  errorBoundMode; //4bits (0.5byte)
    double absErrBound;
    double relBoundRatio;
    double psnr;
    double pw_relBoundRatio;
    int segment_size;
    int pwr_type;
} sz_params;

typedef struct sz_metadata
{
	int versionNumber[3]; //only used for checking the version by calling SZ_GetMetaData()
	int isConstant; //only used for checking if the data are constant values by calling SZ_GetMetaData()
	int isLossless; //only used for checking if the data compression was lossless, used only by calling SZ_GetMetaData()
	int sizeType; //only used for checking whether the size type is "int" or "long" in the compression, used only by calling SZ_GetMetaData()
	size_t dataSeriesLength;
	struct sz_params* conf_params;
} sz_metadata;

extern sz_params *conf_params;

//sz.h
void SZ_Reset();

int SZ_Init(char *configFilePath);

int SZ_Init_Params(sz_params *params);

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

int SZ_compress_args_float_subblock(unsigned char* compressedBytes, float *oriData,
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
size_t s5, size_t s4, size_t s3, size_t s2, size_t s1,
size_t e5, size_t e4, size_t e3, size_t e2, size_t e1,
size_t *outSize, int errBoundMode, double absErr_Bound, double relBoundRatio);

int SZ_compress_args_double_subblock(unsigned char* compressedBytes, double *oriData,
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
size_t s5, size_t s4, size_t s3, size_t s2, size_t s1,
size_t e5, size_t e4, size_t e3, size_t e2, size_t e1,
size_t *outSize, int errBoundMode, double absErr_Bound, double relBoundRatio);

unsigned char *SZ_compress(int dataType, void *data, size_t *outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

unsigned char* SZ_compress_args(int dataType, void *data, size_t *outSize, int errBoundMode, double absErrBound, 
double relBoundRatio, double pwrBoundRatio, int pwrType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

int SZ_compress_args2(int dataType, void *data, unsigned char* compressed_bytes, size_t *outSize, 
int errBoundMode, double absErrBound, double relBoundRatio, double pwrBoundRatio, int pwrType, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

int SZ_compress_args3(int dataType, void *data, unsigned char* compressed_bytes, size_t *outSize, int errBoundMode, double absErrBound, double relBoundRatio, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
size_t s5, size_t s4, size_t s3, size_t s2, size_t s1,
size_t e5, size_t e4, size_t e3, size_t e2, size_t e1);

unsigned char *SZ_compress_rev_args(int dataType, void *data, void *reservedValue, size_t *outSize, int errBoundMode, double absErrBound, double relBoundRatio, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

int SZ_compress_rev_args2(int dataType, void *data, void *reservedValue, unsigned char* compressed_bytes, size_t *outSize, int errBoundMode, double absErrBound, double relBoundRatio, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
unsigned char *SZ_compress_rev(int dataType, void *data, void *reservedValue, size_t *outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void *SZ_decompress(int dataType, unsigned char *bytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
size_t SZ_decompress_args(int dataType, unsigned char *bytes, size_t byteLength, void* decompressed_array, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

sz_metadata* SZ_getMetadata(unsigned char* bytes);
void SZ_printMetadata(sz_metadata* metadata);


void filloutDimArray(size_t* dim, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

size_t compute_total_batch_size();

int isZlibFormat(unsigned char magic1, unsigned char magic2);

unsigned char* SZ_batch_compress(size_t *outSize);

SZ_VarSet* SZ_batch_decompress(unsigned char* compressedStream, size_t length, int *status);

void SZ_Finalize();

void convertSZParamsToBytes(sz_params* params, unsigned char* result);
sz_params* convertBytesToSZParams(unsigned char* bytes);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _SZ_H  ----- */

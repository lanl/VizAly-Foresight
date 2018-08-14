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

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define SZ_VERNUM 0x0130
#define SZ_VER_MAJOR 1
#define SZ_VER_MINOR 4
#define SZ_VER_BUILD 9
#define SZ_VER_REVISION 3

#define HZ 102
#define SZ 101

#define ABS 0
#define REL 1
#define ABS_AND_REL 2
#define ABS_OR_REL 3
#define PW_REL 4

#define SZ_FLOAT 0
#define SZ_DOUBLE 1

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
extern double pw_relBoundRatio;
extern int segment_size;

extern int versionNumber[4];

extern int layers;
extern float predThreshold;
extern int sampleDistance;
extern char optQuantMode;

extern int szMode; //0 (best speed) or 1 (better compression with Gzip)

//extern int spaceFillingCurveTransform; //default is 0, or 1 set by sz.config
//extern int reOrgSize; //the granularity of the reganization of the original data

extern SZ_VarSet* sz_varset;

//typedef unsigned long unsigned long;
//typedef unsigned int uint;

typedef union lshort
{
	unsigned short svalue;
	unsigned char byte[2];
} lshort;

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
	unsigned int max_quant_intervals;
	unsigned int quantization_intervals;
    int dataEndianType;
    int sysEndianType; //sysEndianType is actually set automatically.
    int sol_ID;
    int layers;
    int sampleDistance;
    float predThreshold;    
    int offset;
    int szMode;
    int gzipMode;
    int  errorBoundMode;
    double absErrBound;
    double relBoundRatio;
    double pw_relBoundRatio;
    int segment_size;
} sz_params;

extern sz_params *conf_params;

//conf.c
void updateQuantizationInfo(int quant_intervals);
void clearHuffmanMem();
int SZ_ReadConf();
int SZ_LoadConf();
int checkVersion(char* version);
unsigned int roundUpToPowerOf2(unsigned int base);

//double fabs(double value);

//dataCompression.c
double computeRangeSize_double(double* oriData, int size, double* valueRangeSize, double* medianValue);
float computeRangeSize_float(float* oriData, int size, float* valueRangeSize, float* medianValue);
float computeRangeSize_double_subblock(double* oriData, double* valueRangeSize, double* medianValue,
int r5, int r4, int r3, int r2, int r1,
int s5, int s4, int s3, int s2, int s1,
int e5, int e4, int e3, int e2, int e1);
float computeRangeSize_float_subblock(float* oriData, float* valueRangeSize, float* medianValue,
int r5, int r4, int r3, int r2, int r1,
int s5, int s4, int s3, int s2, int s1,
int e5, int e4, int e3, int e2, int e1);
double min_d(double a, double b);
double max_d(double a, double b);
float min_f(float a, float b);
float max_f(float a, float b);
double getRealPrecision_double(double valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status);
double getRealPrecision_float(float valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status);
void symTransform_8bytes(unsigned char data[8]);
void flush_to_bigEndian_8bytes(unsigned char data[8], int dataEndianType);
void symTransform_2bytes(unsigned char data[2]);
void symTransform_4bytes(unsigned char data[4]);
void flush_to_bigEndian_4bytes(unsigned char data[4]);
void bigEndian_to_OSEndian_double(unsigned char data[8]);
void bigEndian_to_OSEndian_float(unsigned char data[4]);
void compressSingleFloatValue(FloatValueCompressElement *vce, float tgtValue, float precision, float medianValue, 
		int reqLength, int reqBytesLength, int resiBitsLength);
void compressSingleDoubleValue(DoubleValueCompressElement *vce, double tgtValue, double precision, double medianValue, 
		int reqLength, int reqBytesLength, int resiBitsLength);
int compIdenticalLeadingBytesCount_double(unsigned char* preBytes, unsigned char* curBytes);
int compIdenticalLeadingBytesCount_float(unsigned char* preBytes, unsigned char* curBytes);
void addExactData(DynamicByteArray *exactMidByteArray, DynamicIntArray *exactLeadNumArray, 
		DynamicIntArray *resiBitArray, LossyCompressionElement *lce);

//ByteToolkit.c
int bytesToInt_bigEndian(unsigned char* bytes);
void intToBytes_bigEndian(unsigned char *b, unsigned int num);
long bytesToLong_bigEndian(unsigned char* b);
void longToBytes_bigEndian(unsigned char *b, unsigned long num);
long doubleToOSEndianLong(double value);
int floatToOSEndianInt(float value);
short getExponent_float(float value);
short getPrecisionReqLength_float(float precision);
short getExponent_double(double value);
short getPrecisionReqLength_double(double precision);
unsigned char numberOfLeadingZeros_Int(int i);
unsigned char numberOfLeadingZeros_Long(long i);
unsigned char getLeadingNumbers_Int(int v1, int v2);
unsigned char getLeadingNumbers_Long(long v1, long v2);
short bytesToShort(unsigned char* bytes);
int bytesToInt(unsigned char* bytes);
long bytesToLong(unsigned char* bytes);
float bytesToFloat(unsigned char* bytes);
void floatToBytes(unsigned char *b, float num);
double bytesToDouble(unsigned char* bytes);
void doubleToBytes(unsigned char *b, double num);
int extractBytes(unsigned char* byteArray, int k, int validLength);
int getMaskRightCode(int m);
int getLeftMovingCode(int kMod8);
int getRightMovingSteps(int kMod8, int resiBitLength);
int getRightMovingCode(int kMod8, int resiBitLength);
unsigned short* convertByteDataToShortArray(unsigned char* bytes, int byteLength);
void convertShortArrayToBytes(unsigned short* states, int stateLength, unsigned char* bytes);

//TypeManager.c
int convertIntArray2ByteArray_fast_2b(unsigned char* timeStepType, int timeStepTypeLength, unsigned char **result);
void convertByteArray2IntArray_fast_2b(int stepLength, unsigned char* byteArray, int byteArrayLength, unsigned char **intArray);
int convertIntArray2ByteArray_fast_3b(unsigned char* timeStepType, int timeStepTypeLength, unsigned char **result);
void convertByteArray2IntArray_fast_3b(int stepLength, unsigned char* byteArray, int byteArrayLength, unsigned char **intArray);
int getLeftMovingSteps(int k, unsigned char resiBitLength);
int convertIntArray2ByteArray_fast_dynamic(unsigned char* timeStepType, unsigned char* resiBitLength, int resiBitLengthLength, unsigned char **bytes);
int computeBitNumRequired(int dataLength);
void decompressBitArraybySimpleLZ77(int** result, unsigned char* bytes, int bytesLength, int totalLength, int validLength);

//test_zlib.c
unsigned long zlib_compress(unsigned char* data, unsigned long dataLength, unsigned char** compressBytes, int level);
unsigned long zlib_compress2(unsigned char* data, unsigned long dataLength, unsigned char** compressBytes, int level);
unsigned long zlib_compress3(unsigned char* data, unsigned long dataLength, unsigned char* compressBytes, int level);
unsigned long zlib_uncompress(unsigned char* compressBytes, unsigned long cmpSize, unsigned char** oriData, unsigned long targetOriSize);
unsigned long zlib_uncompress2(unsigned char* compressBytes, unsigned long cmpSize, unsigned char** oriData, unsigned long targetOriSize);

//szf.c
void sz_init_c_(char *configFile,int *len,int *ierr);
void sz_finalize_c_();
void SZ_writeData_inBinary_d1_Float_(float* data, char *fileName, int *len);
void sz_compress_d1_float_(float* data, unsigned char *bytes, int *outSize, int *r1);
void sz_compress_d1_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1);
void sz_compress_d2_float_(float* data, unsigned char *bytes, int *outSize, int *r1, int *r2);
void sz_compress_d2_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2);
void sz_compress_d3_float_(float* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3);
void sz_compress_d3_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3);
void sz_compress_d4_float_(float* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d4_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_float_(float* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_compress_d5_float_rev_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_compress_d1_double_(double* data, unsigned char *bytes, int *outSize, int *r1);
void sz_compress_d1_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1);
void sz_compress_d2_double_(double* data, unsigned char *bytes, int *outSize, int *r1, int *r2);
void sz_compress_d2_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2);
void sz_compress_d3_double_(double* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3);
void sz_compress_d3_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3);
void sz_compress_d4_double_(double* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d4_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_double_(double* data, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_compress_d5_double_rev_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_compress_d1_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1);
void sz_compress_d2_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2);
void sz_compress_d3_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3);
void sz_compress_d4_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_float_args_(float* data, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_compress_d1_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1);
void sz_compress_d2_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2);
void sz_compress_d3_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3);
void sz_compress_d4_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_double_args_(double* data, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_compress_d1_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1);
void sz_compress_d2_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2);
void sz_compress_d3_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3);
void sz_compress_d4_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_float_rev_args_(float* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_compress_d1_double_rev_args_(double* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1);
void sz_compress_d2_double_rev_args_(double* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2);
void sz_compress_d3_double_rev_args_(double* data, float *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3);
void sz_compress_d4_double_rev_args_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_compress_d5_double_rev_args_(double* data, double *reservedValue, unsigned char *bytes, int *outSize, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_decompress_d1_float_(unsigned char *bytes, int *byteLength, float *data, int *r1);
void sz_decompress_d2_float_(unsigned char *bytes, int *byteLength, float *data, int *r1, int *r2);
void sz_decompress_d3_float_(unsigned char *bytes, int *byteLength, float *data, int *r1, int *r2, int *r3);
void sz_decompress_d4_float_(unsigned char *bytes, int *byteLength, float *data, int *r1, int *r2, int *r3, int *r4);
void sz_decompress_d5_float_(unsigned char *bytes, int *byteLength, float *data, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_decompress_d1_double_(unsigned char *bytes, int *byteLength, double *data, int *r1);
void sz_decompress_d2_double_(unsigned char *bytes, int *byteLength, double *data, int *r1, int *r2);
void sz_decompress_d3_double_(unsigned char *bytes, int *byteLength, double *data, int *r1, int *r2, int *r3);
void sz_decompress_d4_double_(unsigned char *bytes, int *byteLength, double *data, int *r1, int *r2, int *r3, int *r4);
void sz_decompress_d5_double_(unsigned char *bytes, int *byteLength, double *data, int *r1, int *r2, int *r3, int *r4, int *r5);

void sz_batchaddVar_d1_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1);
void sz_batchaddvar_d2_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2);
void sz_batchaddvar_d3_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3);
void sz_batchaddvar_d4_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_batchaddvar_d5_float_(char* varName, int *len, float* data, int *errBoundMode, float *absErrBound, float *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_batchaddvar_d1_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1);
void sz_batchaddvar_d2_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2);
void sz_batchaddvar_d3_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3);
void sz_batchaddvar_d4_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4);
void sz_batchaddvar_d5_double_(char* varName, int *len, double* data, int *errBoundMode, double *absErrBound, double *relBoundRatio, int *r1, int *r2, int *r3, int *r4, int *r5);
void sz_batchdelvar_c_(char* varName, int *len, int *errState);
void sz_batch_compress_c_(unsigned char* bytes, int *outSize);
void sz_batch_decompress_c_(unsigned char* bytes, int *byteLength, int *ierr);
void sz_getvardim_c_(char* varName, int *len, int *dim, int *r1, int *r2, int *r3, int *r4, int *r5);
void compute_total_batch_size_c_(int *totalSize);
void sz_getvardata_float_(char* varName, int *len, float* data);
void sz_getvardata_double_(char* varName, int *len, double* data);
void sz_freevarset_c_(int *mode);

//sz.h
void SZ_Reset();
int SZ_Init(char *configFilePath);
int SZ_Init_Params(sz_params *params);
int computeDataLength(int r5, int r4, int r3, int r2, int r1);
int computeDimension(int r5, int r4, int r3, int r2, int r1);

void computeReqLength_double(double realPrecision, short radExpo, int* reqLength, double* medianValue);
void computeReqLength_float(float realPrecision, short radExpo, int* reqLength, float* medianValue);

int getPredictionCoefficients(int layers, int dimension, int **coeff_array, int *status);

void compute_segment_precisions_float_1D(float *oriData, int dataLength, float* pwrErrBound, unsigned char* pwrErrBoundBytes);
unsigned int optimize_intervals_float_1D(float *oriData, int dataLength, double realPrecision);
unsigned int optimize_intervals_float_1D_pwr(float *oriData, int dataLength, float* pwrErrBound);

void compute_segment_precisions_float_2D(float *oriData, float* pwrErrBound, 
int r1, int r2, int R2, int edgeSize, unsigned char* pwrErrBoundBytes, float Min, float Max);
unsigned int optimize_intervals_float_2D(float *oriData, int r1, int r2, double realPrecision);
unsigned int optimize_intervals_float_2D_pwr(float *oriData, int r1, int r2, int R2, int edgeSize, float* pwrErrBound);

void compute_segmeoptimize_intervals_double_4Dnt_precisions_float_3D(float *oriData, float* pwrErrBound,
int r1, int r2, int r3, int R2, int R3, int edgeSize, unsigned char* pwrErrBoundBytes, float Min, float Max);
unsigned int optimize_intervals_float_3D(float *oriData, int r1, int r2, int r3, double realPrecision);
unsigned int optimize_intervals_float_3D_pwr(float *oriData, int r1, int r2, int r3, int R2, int R3, int edgeSize, float* pwrErrBound);

unsigned int optimize_intervals_float_4D(float *oriData, int r1, int r2, int r3, int r4, double realPrecision);

void compute_segment_precisions_double_1D(double *oriData, int dataLength, double* pwrErrBound, unsigned char* pwrErrBoundBytes);
unsigned int optimize_intervals_double_1D_pwr(double *oriData, int dataLength, double* pwrErrBound);
unsigned int optimize_intervals_double_1D(double *oriData, int dataLength, double realPrecision);

void compute_segment_precisions_double_2D(double *oriData, double* pwrErrBound, 
int r1, int r2, int R2, int edgeSize, unsigned char* pwrErrBoundBytes, double Min, double Max);
unsigned int optimize_intervals_double_2D(double *oriData, int r1, int r2, double realPrecision);
unsigned int optimize_intervals_double_2D_pwr(double *oriData, int r1, int r2, int R2, int edgeSize, double* pwrErrBound);

void compute_segment_precisions_double_3D(double *oriData, double* pwrErrBound, 
int r1, int r2, int r3, int R2, int R3, int edgeSize, unsigned char* pwrErrBoundBytes, double Min, double Max);
unsigned int optimize_intervals_double_3D(double *oriData, int r1, int r2, int r3, double realPrecision);
unsigned int optimize_intervals_double_3D_pwr(double *oriData, int r1, int r2, int r3, int R2, int R3, int edgeSize, double* pwrErrBound);

unsigned int optimize_intervals_double_4D(double *oriData, int r1, int r2, int r3, int r4, double realPrecision);

unsigned int optimize_intervals_float_1D_subblock(float *oriData, double realPrecision, int r1, int s1, int e1);
unsigned int optimize_intervals_float_2D_subblock(float *oriData, double realPrecision, int r1, int r2, int s1, int s2, int e1, int e2);
unsigned int optimize_intervals_float_3D_subblock(float *oriData, double realPrecision, int r1, int r2, int r3, int s1, int s2, int s3, int e1, int e2, int e3);
unsigned int optimize_intervals_float_4D_subblock(float *oriData, double realPrecision, int r1, int r2, int r3, int r4, int s1, int s2, int s3, int s4, int e1, int e2, int e3, int e4);


unsigned int optimize_intervals_double_1D_subblock(double *oriData, double realPrecision, int r1, int s1, int e1);;
unsigned int optimize_intervals_double_2D_subblock(double *oriData, double realPrecision, int r1, int r2, int s1, int s2, int e1, int e2);
unsigned int optimize_intervals_double_3D_subblock(double *oriData, double realPrecision, int r1, int r2, int r3, int s1, int s2, int s3, int e1, int e2, int e3);
unsigned int optimize_intervals_double_4D_subblock(double *oriData, double realPrecision, int r1, int r2, int r3, int r4, int s1, int s2, int s3, int s4, int e1, int e2, int e3, int e4);

int computeBlockEdgeSize(int segmentSize);
int computeBlockEdgeSize_3D(int segmentSize);
int computeBlockEdgeSize_2D(int segmentSize);

void SZ_compress_args_float_StoreOriData(float* oriData, int dataLength, TightDataPointStorageF* tdps,
unsigned char** newByteData, int *outSize);

TightDataPointStorageF* SZ_compress_float_1D_MDQ(float *oriData, 
int dataLength, double realPrecision, float valueRangeSize, float medianValue_f);
void SZ_compress_args_float_NoCkRngeNoGzip_1D(unsigned char** newByteData, float *oriData, int dataLength, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f);
void SZ_compress_args_float_NoCkRngeNoGzip_1D_pwr(unsigned char** newByteData, float *oriData, int dataLength, 
int *outSize, float min, float max);

TightDataPointStorageF* SZ_compress_float_2D_MDQ(float *oriData, int r1, int r2, double realPrecision, float valueRangeSize, float medianValue_f);
void SZ_compress_args_float_NoCkRngeNoGzip_2D(unsigned char** newByteData, float *oriData, int r1, int r2, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f);
void SZ_compress_args_float_NoCkRngeNoGzip_2D_pwr(unsigned char** newByteData, float *oriData, int r1, int r2, 
int *outSize, float min, float max);

TightDataPointStorageF* SZ_compress_float_3D_MDQ(float *oriData, int r1, int r2, int r3, double realPrecision, float valueRangeSize, float medianValue_f);
void SZ_compress_args_float_NoCkRngeNoGzip_3D(unsigned char** newByteData, float *oriData, int r1, int r2, int r3, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f);
void SZ_compress_args_float_NoCkRngeNoGzip_3D_pwr(unsigned char** newByteData, float *oriData, int r1, int r2, int r3, 
int *outSize, float min, float max);

TightDataPointStorageF* SZ_compress_float_4D_MDQ(float *oriData, int r1, int r2, int r3, int r4, double realPrecision, float valueRangeSize, float medianValue_f);
void SZ_compress_args_float_NoCkRngeNoGzip_4D(unsigned char** newByteData, float *oriData, int r1, int r2, int r3, int r4, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f);


TightDataPointStorageD* SZ_compress_double_1D_MDQ(double *oriData, 
int dataLength, double realPrecision, double valueRangeSize, double medianValue_d);
void SZ_compress_args_double_StoreOriData(double* oriData, int dataLength, TightDataPointStorageD* tdps, 
unsigned char** newByteData, int *outSize);
void SZ_compress_args_double_NoCkRngeNoGzip_1D(unsigned char** newByteData, double *oriData, int dataLength, double realPrecision, int *outSize, double valueRangeSize, double medianValue_d);
void SZ_compress_args_double_NoCkRngeNoGzip_1D_pwr(unsigned char** newByteData, double *oriData, 
int dataLength, int *outSize, double min, double max);

TightDataPointStorageD* SZ_compress_double_2D_MDQ(double *oriData, int r1, int r2, double realPrecision, double valueRangeSize, double medianValue_d);
void SZ_compress_args_double_NoCkRngeNoGzip_2D(unsigned char** newByteData, double *oriData, int r1, int r2, double realPrecision, int *outSize, double valueRangeSize, double medianValue_d);
void SZ_compress_args_double_NoCkRngeNoGzip_2D_pwr(unsigned char** newByteData, double *oriData, int r1, int r2, int *outSize, double min, double max);

TightDataPointStorageD* SZ_compress_double_3D_MDQ(double *oriData, int r1, int r2, int r3, double realPrecision, double valueRangeSize, double medianValue_d);
void SZ_compress_args_double_NoCkRngeNoGzip_3D(unsigned char** newByteData, double *oriData, int r1, int r2, int r3, double realPrecision, int *outSize, double valueRangeSize, double medianValue_d);
void SZ_compress_args_double_NoCkRngeNoGzip_3D_pwr(unsigned char** newByteData, double *oriData, int r1, int r2, int r3, 
int *outSize, double min, double max);

TightDataPointStorageD* SZ_compress_double_4D_MDQ(double *oriData, int r1, int r2, int r3, int r4, double realPrecision, double valueRangeSize, double medianValue_d);
void SZ_compress_args_double_NoCkRngeNoGzip_4D(unsigned char** newByteData, double *oriData, int r1, int r2, int r3, int r4, double realPrecision, int *outSize, double valueRangeSize, double medianValue_d);


void SZ_compress_args_float_NoCkRnge_1D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f,
int r1, int s1, int e1);
void SZ_compress_args_float_NoCkRnge_2D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f,
int r2, int r1, int s2, int s1, int e2, int e1);
void SZ_compress_args_float_NoCkRnge_3D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f,
int r3, int r2, int r1, int s3, int s2, int s1, int e3, int e2, int e1);
void SZ_compress_args_float_NoCkRnge_4D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, int *outSize, float valueRangeSize, float medianValue_f,
int r4, int r3, int r2, int r1, int s4, int s3, int s2, int s1, int e4, int e3, int e2, int e1);

void SZ_compress_args_double_NoCkRnge_1D_subblock(unsigned char* compressedBytes, double *oriData, double realPrecision, int *outSize, double valueRangeSize, double medianValue_d,
int r1, int s1, int e1);
void SZ_compress_args_double_NoCkRnge_2D_subblock(unsigned char* compressedBytes, double *oriData, double realPrecision, int *outSize, double valueRangeSize, double medianValue_d,
int r2, int r1, int s2, int s1, int e2, int e1);
void SZ_compress_args_double_NoCkRnge_3D_subblock(unsigned char* compressedBytes, double *oriData, double realPrecision, int *outSize, double valueRangeSize, double medianValue_d,
int r3, int r2, int r1, int s3, int s2, int s1, int e3, int e2, int e1);
void SZ_compress_args_double_NoCkRnge_4D_subblock(unsigned char* compressedBytes, double *oriData, double realPrecision, int *outSize, double valueRangeSize, double medianValue_d,
int r4, int r3, int r2, int r1, int s4, int s3, int s2, int s1, int e4, int e3, int e2, int e1);


TightDataPointStorageF* SZ_compress_float_1D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
int r1, int s1, int e1);
TightDataPointStorageF* SZ_compress_float_2D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
int r1, int r2, int s1, int s2, int e1, int e2);
TightDataPointStorageF* SZ_compress_float_3D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
int r1, int r2, int r3, int s1, int s2, int s3, int e1, int e2, int e3);
TightDataPointStorageF* SZ_compress_float_4D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
int r1, int r2, int r3, int r4, int s1, int s2, int s3, int s4, int e1, int e2, int e3, int e4);

TightDataPointStorageD* SZ_compress_double_1D_MDQ_subblock(double *oriData, double realPrecision, double valueRangeSize, double medianValue_d,
int r1, int s1, int e1);
TightDataPointStorageD* SZ_compress_double_2D_MDQ_subblock(double *oriData, double realPrecision, double valueRangeSize, double medianValue_d,
int r1, int r2, int s1, int s2, int e1, int e2);
TightDataPointStorageD* SZ_compress_double_3D_MDQ_subblock(double *oriData, double realPrecision, double valueRangeSize, double medianValue_d,
int r1, int r2, int r3, int s1, int s2, int s3, int e1, int e2, int e3);
TightDataPointStorageD* SZ_compress_double_4D_MDQ_subblock(double *oriData, double realPrecision, double valueRangeSize, double medianValue_d,
int r1, int r2, int r3, int r4, int s1, int s2, int s3, int s4, int e1, int e2, int e3, int e4);

void SZ_compress_args_float_withinRange(unsigned char** newByteData, float *oriData, int dataLength, int *outSize);
void SZ_compress_args_double_withinRange(unsigned char** newByteData, double *oriData, int dataLength, int *outSize);

int SZ_compress_args_float(unsigned char** newByteData, float *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio);
int SZ_compress_args_double(unsigned char** newByteData, double *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio);

int SZ_compress_args_float_subblock(unsigned char* compressedBytes, float *oriData, 
int r5, int r4, int r3, int r2, int r1,
int s5, int s4, int s3, int s2, int s1,
int e5, int e4, int e3, int e2, int e1,
int *outSize, int errBoundMode, double absErr_Bound, double relBoundRatio);

int SZ_compress_args_double_subblock(unsigned char* compressedBytes, double *oriData, 
int r5, int r4, int r3, int r2, int r1,
int s5, int s4, int s3, int s2, int s1,
int e5, int e4, int e3, int e2, int e1,
int *outSize, int errBoundMode, double absErr_Bound, double relBoundRatio);

int SZ_compress_args_float_wRngeNoGzip(unsigned char** newByteData, float *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio);
int SZ_compress_args_double_wRngeNoGzip(unsigned char** newByteData, double *oriData, 
int r5, int r4, int r3, int r2, int r1, int *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio);

unsigned char *SZ_compress(int dataType, void *data, int *outSize, int r5, int r4, int r3, int r2, int r1);
unsigned char *SZ_compress_args(int dataType, void *data, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1);
int SZ_compress_args2(int dataType, void *data, unsigned char* compressed_bytes, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1);
int SZ_compress_args3(int dataType, void *data, unsigned char* compressed_bytes, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, 
int r5, int r4, int r3, int r2, int r1, 
int s5, int s4, int s3, int s2, int s1,
int e5, int e4, int e3, int e2, int e1);

unsigned char *SZ_compress_rev_args(int dataType, void *data, void *reservedValue, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1);
int SZ_compress_rev_args2(int dataType, void *data, void *reservedValue, unsigned char* compressed_bytes, int *outSize, int errBoundMode, double absErrBound, double relBoundRatio, int r5, int r4, int r3, int r2, int r1);
unsigned char *SZ_compress_rev(int dataType, void *data, void *reservedValue, int *outSize, int r5, int r4, int r3, int r2, int r1);

int SZ_decompress_args_float(float** newData, int r5, int r4, int r3, int r2, int r1, unsigned char* cmpBytes, int cmpSize);
int SZ_decompress_args_double(double** newData, int r5, int r4, int r3, int r2, int r1, unsigned char* cmpBytes, int cmpSize);
void *SZ_decompress(int dataType, unsigned char *bytes, int byteLength, int r5, int r4, int r3, int r2, int r1);
int SZ_decompress_args(int dataType, unsigned char *bytes, int byteLength, void* decompressed_array, int r5, int r4, int r3, int r2, int r1);

void filloutDimArray(int* dim, int r5, int r4, int r3, int r2, int r1);
int compute_total_batch_size();
int isZlibFormat(unsigned char magic1, unsigned char magic2);

unsigned char* SZ_batch_compress(int *outSize);
SZ_VarSet* SZ_batch_decompress(unsigned char* compressedStream, int length, int *status);

void SZ_Finalize();

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _SZ_H  ----- */

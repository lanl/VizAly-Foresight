/**
 *  @file sz_float.h
 *  @author Sheng Di
 *  @date July, 2017
 *  @brief Header file for the sz_float.c.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _SZ_Float_H
#define _SZ_Float_H

#ifdef __cplusplus
extern "C" {
#endif

void computeReqLength_float(float realPrecision, short radExpo, int* reqLength, float* medianValue);

unsigned int optimize_intervals_float_1D(float *oriData, size_t dataLength, double realPrecision);
unsigned int optimize_intervals_float_2D(float *oriData, size_t r1, size_t r2, double realPrecision);
unsigned int optimize_intervals_float_3D(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision);
unsigned int optimize_intervals_float_4D(float *oriData, size_t r1, size_t r2, size_t r3, size_t r4, double realPrecision);

TightDataPointStorageF* SZ_compress_float_1D_MDQ(float *oriData, 
size_t dataLength, double realPrecision, float valueRangeSize, float medianValue_f);

void SZ_compress_args_float_StoreOriData(float* oriData, size_t dataLength, TightDataPointStorageF* tdps, 
unsigned char** newByteData, size_t *outSize);

void SZ_compress_args_float_NoCkRngeNoGzip_1D(unsigned char** newByteData, float *oriData, 
size_t dataLength, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f);

TightDataPointStorageF* SZ_compress_float_2D_MDQ(float *oriData, size_t r1, size_t r2, double realPrecision, float valueRangeSize, float medianValue_f);

void SZ_compress_args_float_NoCkRngeNoGzip_2D(unsigned char** newByteData, float *oriData, size_t r1, size_t r2, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f);

TightDataPointStorageF* SZ_compress_float_3D_MDQ(float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, float valueRangeSize, float medianValue_f);

void SZ_compress_args_float_NoCkRngeNoGzip_3D(unsigned char** newByteData, float *oriData, size_t r1, size_t r2, size_t r3, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f);

TightDataPointStorageF* SZ_compress_float_4D_MDQ(float *oriData, size_t r1, size_t r2, size_t r3, size_t r4, double realPrecision, float valueRangeSize, float medianValue_f);

void SZ_compress_args_float_NoCkRngeNoGzip_4D(unsigned char** newByteData, float *oriData, size_t r1, size_t r2, size_t r3, size_t r4, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f);

void SZ_compress_args_float_withinRange(unsigned char** newByteData, float *oriData, size_t dataLength, size_t *outSize);

int SZ_compress_args_float_wRngeNoGzip(unsigned char** newByteData, float *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio, double pwrErrRatio);

int SZ_compress_args_float(unsigned char** newByteData, float *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio, double pwRelBoundRatio, int pwrType);

int SZ_compress_args_float_subblock(unsigned char* compressedBytes, float *oriData,
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
size_t s5, size_t s4, size_t s3, size_t s2, size_t s1,
size_t e5, size_t e4, size_t e3, size_t e2, size_t e1,
size_t *outSize, int errBoundMode, double absErr_Bound, double relBoundRatio);

void SZ_compress_args_float_NoCkRnge_1D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f,
size_t r1, size_t s1, size_t e1); 

void SZ_compress_args_float_NoCkRnge_2D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f,
size_t r2, size_t r1, size_t s2, size_t s1, size_t e2, size_t e1); 

void SZ_compress_args_float_NoCkRnge_3D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f,
size_t r3, size_t r2, size_t r1, size_t s3, size_t s2, size_t s1, size_t e3, size_t e2, size_t e1); 

void SZ_compress_args_float_NoCkRnge_4D_subblock(unsigned char* compressedBytes, float *oriData, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f,
size_t r4, size_t r3, size_t r2, size_t r1, size_t s4, size_t s3, size_t s2, size_t s1, size_t e4, size_t e3, size_t e2, size_t e1);

unsigned int optimize_intervals_float_1D_subblock(float *oriData, double realPrecision, size_t r1, size_t s1, size_t e1); 
unsigned int optimize_intervals_float_2D_subblock(float *oriData, double realPrecision, size_t r1, size_t r2, size_t s1, size_t s2, size_t e1, size_t e2); 
unsigned int optimize_intervals_float_3D_subblock(float *oriData, double realPrecision, size_t r1, size_t r2, size_t r3, size_t s1, size_t s2, size_t s3, size_t e1, size_t e2, size_t e3); 
unsigned int optimize_intervals_float_4D_subblock(float *oriData, double realPrecision, size_t r1, size_t r2, size_t r3, size_t r4, size_t s1, size_t s2, size_t s3, size_t s4, size_t e1, size_t e2, size_t e3, size_t e4);

TightDataPointStorageF* SZ_compress_float_1D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
size_t r1, size_t s1, size_t e1); 

TightDataPointStorageF* SZ_compress_float_2D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
size_t r1, size_t r2, size_t s1, size_t s2, size_t e1, size_t e2); 

TightDataPointStorageF* SZ_compress_float_3D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
size_t r1, size_t r2, size_t r3, size_t s1, size_t s2, size_t s3, size_t e1, size_t e2, size_t e3); 

TightDataPointStorageF* SZ_compress_float_4D_MDQ_subblock(float *oriData, double realPrecision, float valueRangeSize, float medianValue_f,
size_t r1, size_t r2, size_t r3, size_t r4, size_t s1, size_t s2, size_t s3, size_t s4, size_t e1, size_t e2, size_t e3, size_t e4);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _SZ_Float_H  ----- */


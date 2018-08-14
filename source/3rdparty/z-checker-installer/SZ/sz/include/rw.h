/**
 *  @file io.h
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief Header file for the whole io interface.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _IO_H
#define _IO_H

#include <stdio.h>

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

float** create2DArray_float(int m, int n);
float*** create3DArray_float(int p, int m, int n);
void free2DArray_float(float** data, int m);
void free3DArray_float(float*** data, int p, int m);

double** create2DArray_double(int m, int n);
double*** create3DArray_double(int p, int m, int n);
void free2DArray_double(double** data, int m);
void free3DArray_double(double*** data, int p, int m);

int checkFileSize(char *srcFilePath, int *status);
unsigned char *readByteData(char *srcFilePath, int *byteLength, int *status);
double *readDoubleData_systemEndian(char *srcFilePath, int *nbEle, int *status);
float *readFloatData_systemEndian(char *srcFilePath, int *nbEle, int *status);
double *readDoubleData(char *srcFilePath, int *nbEle, int *status);
float *readFloatData(char *srcFilePath, int *nbEle, int *status);
void writeByteData(unsigned char *bytes, int outSize, char *tgtFilePath, int *status);
void writeDoubleData(double *data, int nbEle, char *tgtFilePath, int *status);
void writeFloatData(float *data, int nbEle, char *tgtFilePath, int *status);
void writeData(void *data, int dataType, int nbEle, char *tgtFilePath, int *status);
void writeFloatData_inBytes(float *data, int nbEle, char* tgtFilePath, int *status);
void writeDoubleData_inBytes(double *data, int nbEle, char* tgtFilePath, int *status);
void writeShortData(unsigned short *states, int stateLength, char *tgtFilePath, int *status);
unsigned short* readShortData(char *srcFilePath, int *dataLength, int *status);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _IO_H  ----- */

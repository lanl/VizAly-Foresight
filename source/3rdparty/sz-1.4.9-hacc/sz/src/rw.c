/**
 *  @file rw.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief io interface for fortrance
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rw.h"
#include "sz.h"

float** create2DArray_float(int m, int n)
{
	int i=0;
	float **data = (float**)malloc(sizeof(float*)*m);
	for(i=0;i<m;i++)
		data[i] = (float*)malloc(sizeof(float)*n);
	return data;
}

void free2DArray_float(float** data, int m)
{
	int i = 0;
	for(i=0;i<m;i++)
		free(data[i]);
	free(data);	
}

float*** create3DArray_float(int p, int m, int n)
{
	int i = 0, j = 0;
	float ***data = (float***)malloc(sizeof(float**)*m);
	for(i=0;i<p;i++)
	{
		data[i] = (float**)malloc(sizeof(float*)*n);
		for(j=0;j<m;j++)
			data[i][j] = (float*)malloc(sizeof(float)*n);
	}
	return data;
}

void free3DArray_float(float*** data, int p, int m)
{
	int i,j;
	for(i=0;i<p;i++)
	{
		for(j=0;j<m;j++)
			free(data[i][j]);
		free(data[i]);
	}
	free(data);	
}

double** create2DArray_double(int m, int n)
{
	int i=0;
	double **data = (double**)malloc(sizeof(double*)*m);
	for(i=0;i<m;i++)
			data[i] = (double*)malloc(sizeof(double)*n);
			
	return data;
}

void free2DArray_double(double** data, int m)
{
	int i;
	for(i=0;i<m;i++)
		free(data[i]);
	free(data);	
}

double*** create3DArray_double(int p, int m, int n)
{
	int i = 0, j = 0;
	double ***data = (double***)malloc(sizeof(double**)*m);
	for(i=0;i<p;i++)
	{
		data[i] = (double**)malloc(sizeof(double*)*n);
		for(j=0;j<m;j++)
			data[i][j] = (double*)malloc(sizeof(double)*n);
	}
	return data;
}

void free3DArray_double(double*** data, int p, int m)
{
	int i,j;
	for(i=0;i<p;i++)
	{
		for(j=0;j<m;j++)
			free(data[i][j]);
		free(data[i]);
	}
	free(data);	
}

int checkFileSize(char *srcFilePath, int *status)
{
	int filesize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
	{
		printf("Failed to open input file. 1\n");
		*status = SZ_FERR;
		return -1;
	}
	fseek(pFile, 0, SEEK_END);
    filesize = (int)ftell(pFile);
    fclose(pFile);
    *status = SZ_SCES;
    return filesize;
}

unsigned char *readByteData(char *srcFilePath, int *byteLength, int *status)
{
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        *status = SZ_FERR;
        return 0;
    }
	fseek(pFile, 0, SEEK_END);
    *byteLength = (int)ftell(pFile);
    fclose(pFile);
    
    unsigned char *byteBuf = ( unsigned char *)malloc((*byteLength)*sizeof(unsigned char)); //sizeof(char)==1
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        *status = SZ_FERR;
        return 0;
    }
    fread(byteBuf, 1, *byteLength, pFile);
    fclose(pFile);
    *status = SZ_SCES;
    return byteBuf;
}

double *readDoubleData(char *srcFilePath, int *nbEle, int *status)
{
	int state = SZ_SCES;
	if(dataEndianType==sysEndianType)
	{
		double *daBuf = readDoubleData_systemEndian(srcFilePath, nbEle,&state);
		*status = state;
		return daBuf;
	}
	else
	{
		int i,j;
		
		int byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state==SZ_FERR)
		{
			*status = SZ_FERR;
			return NULL;
		}
		double *daBuf = (double *)malloc(byteLength);
		*nbEle = byteLength/8;
		
		ldouble buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*8;
			memcpy(buf.byte, bytes+j, 8);
			symTransform_8bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		return daBuf;
	}
}

float *readFloatData(char *srcFilePath, int *nbEle, int *status)
{
	int state = SZ_SCES;
	if(dataEndianType==sysEndianType)
	{
		float *daBuf = readFloatData_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		int i,j;
		
		int byteLength;
		unsigned char* bytes = readByteData(srcFilePath, &byteLength, &state);
		if(state == SZ_FERR)
		{
			*status = SZ_FERR;
			return NULL;
		}
		float *daBuf = (float *)malloc(byteLength);
		*nbEle = byteLength/4;
		
		lfloat buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*4;
			memcpy(buf.byte, bytes+j, 4);
			symTransform_4bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		return daBuf;
	}
}

double *readDoubleData_systemEndian(char *srcFilePath, int *nbEle, int *status)
{
	int inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        *status = SZ_FERR;
        return NULL;
    }
	fseek(pFile, 0, SEEK_END);
    inSize = ftell(pFile);
    *nbEle = (int)inSize/8; //only support double in this version
    fclose(pFile);
    
    double *daBuf = (double *)malloc(inSize);
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        *status = SZ_FERR;
        return NULL;
    }
    fread(daBuf, 8, *nbEle, pFile);
    fclose(pFile);
    *status = SZ_SCES;
    return daBuf;
}

float *readFloatData_systemEndian(char *srcFilePath, int *nbEle, int *status)
{
	int inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        *status = SZ_FERR;
        return NULL;
    }
	fseek(pFile, 0, SEEK_END);
    inSize = ftell(pFile);
    *nbEle = (int)inSize/4; 
    fclose(pFile);
    
    if(inSize<=0)
    {
		printf("Error: input file is wrong!\n");
		*status = SZ_FERR;
	}
    
    float *daBuf = (float *)malloc(inSize);
    
    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        *status = SZ_FERR;
        return NULL;
    }
    fread(daBuf, 4, *nbEle, pFile);
    fclose(pFile);
    *status = SZ_SCES;
    return daBuf;
}

void writeByteData(unsigned char *bytes, int byteLength, char *tgtFilePath, int *status)
{
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        *status = SZ_FERR;
        return;
    }
    
    fwrite(bytes, 1, byteLength, pFile); //write outSize bytes
    fclose(pFile);
    *status = SZ_SCES;
}

void writeDoubleData(double *data, int nbEle, char *tgtFilePath, int *status)
{
	int i = 0;
	char s[64];
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        *status = SZ_FERR;
        return;
    }
    
    for(i = 0;i<nbEle;i++)
	{
		sprintf(s,"%.20G\n",data[i]);
		fputs(s, pFile);
	}
    
    fclose(pFile);
    *status = SZ_SCES;
}

void writeFloatData(float *data, int nbEle, char *tgtFilePath, int *status)
{
	int i = 0;
	char s[64];
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        *status = SZ_FERR;
        return;
    }
   
    for(i = 0;i<nbEle;i++)
	{
		//printf("i=%d\n",i);
		//printf("data[i]=%f\n",data[i]);
		sprintf(s,"%.30G\n",data[i]);
		fputs(s, pFile);
	}
    
    fclose(pFile);
    *status = SZ_SCES;
}

void writeData(void *data, int dataType, int nbEle, char *tgtFilePath, int *status)
{
	int state = SZ_SCES;
	if(dataType == SZ_FLOAT)
	{
		float* dataArray = (float *)data;
		writeFloatData(dataArray, nbEle, tgtFilePath, &state);
	}
	else if(dataType == SZ_DOUBLE)
	{
		double* dataArray = (double *)data;
		writeDoubleData(dataArray, nbEle, tgtFilePath, &state);	
	}
	else
	{
		printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
		*status = SZ_TERR; //wrong type
		return;
	}
	*status = state;
}

void writeFloatData_inBytes(float *data, int nbEle, char* tgtFilePath, int *status)
{
	int i = 0, state = SZ_SCES;
	lfloat buf;
	unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(float));
	for(i=0;i<nbEle;i++)
	{
		buf.value = data[i];
		bytes[i*4+0] = buf.byte[0];
		bytes[i*4+1] = buf.byte[1];
		bytes[i*4+2] = buf.byte[2];
		bytes[i*4+3] = buf.byte[3];					
	}

	int byteLength = nbEle*sizeof(float);
	writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void writeDoubleData_inBytes(double *data, int nbEle, char* tgtFilePath, int *status)
{
	int i = 0, index = 0, state = SZ_SCES;
	ldouble buf;
	unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(double));
	for(i=0;i<nbEle;i++)
	{
		index = i*8;
		buf.value = data[i];
		bytes[index+0] = buf.byte[0];
		bytes[index+1] = buf.byte[1];
		bytes[index+2] = buf.byte[2];
		bytes[index+3] = buf.byte[3];
		bytes[index+4] = buf.byte[4];
		bytes[index+5] = buf.byte[5];
		bytes[index+6] = buf.byte[6];
		bytes[index+7] = buf.byte[7];
	}

	int byteLength = nbEle*sizeof(double);
	writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void writeShortData(unsigned short *states, int stateLength, char *tgtFilePath, int *status)
{
	int state = SZ_SCES;
	int byteLength = stateLength*2;
	unsigned char* bytes = (unsigned char*)malloc(byteLength*sizeof(char));
	convertShortArrayToBytes(states, stateLength, bytes);
	writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

unsigned short* readShortData(char *srcFilePath, int *dataLength, int *status)
{
	int byteLength = 0, state = SZ_SCES;
	unsigned char * bytes = readByteData(srcFilePath, &byteLength, &state);
	*dataLength = byteLength/2;
	unsigned short* states = convertByteDataToShortArray(bytes, byteLength);
	free(bytes);
	*status = state;
	return states;
}

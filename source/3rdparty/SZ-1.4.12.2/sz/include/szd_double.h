/**
 *  @file szd_double.h
 *  @author Sheng Di
 *  @date July, 2017
 *  @brief Header file for the szd_double.c.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _SZD_Double_H
#define _SZD_Double_H

#ifdef __cplusplus
extern "C" {
#endif

#include "TightDataPointStorageD.h"

void decompressDataSeries_double_1D(double** data, size_t dataSeriesLength, TightDataPointStorageD* tdps);
void decompressDataSeries_double_2D(double** data, size_t r1, size_t r2, TightDataPointStorageD* tdps);
void decompressDataSeries_double_3D(double** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageD* tdps);
void decompressDataSeries_double_4D(double** data, size_t r1, size_t r2, size_t r3, size_t r4, TightDataPointStorageD* tdps);
void getSnapshotData_double_1D(double** data, size_t dataSeriesLength, TightDataPointStorageD* tdps, int errBoundMode);
void getSnapshotData_double_2D(double** data, size_t r1, size_t r2, TightDataPointStorageD* tdps, int errBoundMode);
void getSnapshotData_double_3D(double** data, size_t r1, size_t r2, size_t r3, TightDataPointStorageD* tdps, int errBoundMode);
void getSnapshotData_double_4D(double** data, size_t r1, size_t r2, size_t r3, size_t r4, TightDataPointStorageD* tdps, int errBoundMode);

int SZ_decompress_args_double(double** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _SZD_Double_H  ----- */

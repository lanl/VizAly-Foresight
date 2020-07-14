/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/
#ifdef CBENCH_HAS_ZFP

#ifndef _ZFP_COMPRESSOR_H_
#define _ZFP_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
#include "fftw3.h"

#include <zfp.h>

#define REAL 0
#define IMAG 1
#define PI 3.14159265358979

class ZFPCompressor: public CompressorInterface
{
    size_t zfpCompressedSize;
    int numDims;

  public:
    ZFPCompressor();
    ~ZFPCompressor();

    void init();
    int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    void close();

    zfp_type getZfpType(std::string dataType);
};


inline ZFPCompressor::ZFPCompressor()
{
    compressorName = "zfp";
    numDims = 1;
}

inline ZFPCompressor::~ZFPCompressor()
{
    
}


inline void ZFPCompressor::init()
{

}


inline zfp_type ZFPCompressor::getZfpType(std::string dataType)
{
	if (dataType == "float")
		return zfp_type_float;
	else if (dataType == "double")
		return zfp_type_double;
	else if ((dataType == "int32_t") || (dataType == "int16_t"))
		return zfp_type_int32;
	else if (dataType == "int64_t")
		return zfp_type_int64;
}


inline int ZFPCompressor::compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	size_t numel = n[0];
    int numDims = 1;
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
        {
			numel *= n[i];
            numDims++;
        }

    // Read in json compression parameters
    double abs = 1E-3;
    int rel = 32;
    bool compressionTypeAbs = true;
    if ( compressorParameters.find("abs") != compressorParameters.end() )
    {
        abs = strConvert::to_double( compressorParameters["abs"] );
        compressionTypeAbs = true;
    }
    else if ( compressorParameters.find("rel") != compressorParameters.end())
    {
        rel = strConvert::to_int( compressorParameters["rel"] );
        compressionTypeAbs = false;
    }




        Timer entropyClock;
        entropyClock.start();
        float entropy = 0;
        float value_max = 0;
        float average = 0;
        float value_min = 3.402823466e+38F;
        int pos[32];
        for (int i = 0; i < 32; i++)
                pos[i] = 0;
        for (int i = 0; i < numel; i++) {
                average += static_cast<float *>(input)[i];
                if (static_cast<float *>(input)[i] > value_max)
                        value_max = static_cast<float *>(input)[i];
                if (static_cast<float *>(input)[i] < value_min)
                        value_min = static_cast<float *>(input)[i];
        }
        float minmax = value_max - value_min;
//        printf("min = %f, max = %f, minmax = %f, numel = %d\n", value_min, value_max, minmax, numel);
        for (int i = 0; i < numel; i++)
                pos[(int)floor(((static_cast<float *>(input)[i] - value_min)/minmax) * 32)] += 1;
        for (int i = 0; i < 32; i++) {
                if (pos[i] != 0)
                        entropy += (-float(pos[i])/numel) * (std::log(float(pos[i])/numel)/std::log(2));
        }
        entropyClock.stop();
        double entropy_time = entropyClock.getDuration();
        average = average / numel;
        printf("%f ", entropy_time);
        printf("%f ", entropy);
        printf("%f ", average);

/*
        Timer lorenzoClock;
        lorenzoClock.start();
        float lorenzo = 0;
        int helper0000 = 0;
        float msr_lorenzo = 0;
        float layer01, layer02, layer03;
        for (int i = n[0]*n[0]+n[0]; i < numel; i++) {
                if (i%(n[0]*n[0]) == 0)
                        i += n[0];
                if (i%n[0] == 0)
                        i += 1;
                helper0000 += 1;
                layer03 = static_cast<float *>(input)[i-n[0]*n[1]-n[0]-1];
                layer02 = -static_cast<float *>(input)[i-n[0]*n[1]-n[0]] - static_cast<float *>(input)[i-n[0]*n[1]-1] - static_cast<float *>(input)[i-n[0]-1];
                layer01 = static_cast<float *>(input)[i-n[0]*n[1]] + static_cast<float *>(input)[i-n[0]] + static_cast<float *>(input)[i-1];
                msr_lorenzo += pow(layer01+layer02+layer03-static_cast<float *>(input)[i], 2);
                if (fabs(layer01+layer02+layer03-static_cast<float *>(input)[i]) > absTol*127) {
//                      if (layer01+layer02+layer03-static_cast<float *>(input)[i] > 0)
//                              static_cast<float *>(input)[i] = layer01+layer02+layer03-(absTol*127);
//                      else
//                              static_cast<float *>(input)[i] = layer01+layer02+layer03+(absTol*127);
                        lorenzo += 1;
                }
        }
//      printf("%d ", helper0000);
        msr_lorenzo = sqrt(msr_lorenzo/helper0000);
        lorenzo = lorenzo/helper0000;
        lorenzoClock.stop();
        double lorenzo_time = lorenzoClock.getDuration();
        printf("%f ", lorenzo_time);
        printf("%f ", lorenzo);
        printf("%f ", msr_lorenzo);
*/


	Timer fft3dClock;
	fft3dClock.start();
        int len = n[0];
        int lena = len/2+1;
        fftwf_complex y[len*len*len];
        float z[len*len*len];
	float energy_low, energy_mid_low, energy_mid_high, energy_high;
/*
        for (int i = 0; i < len; i++) {
                for (int j = 0; j < len; j++) {
                        for (int k = 0; k < len; k++) {
                                x[i*len*len+j*len+k] = sin((-i+j+k)*2*PI/len*2+22);
                //x[i] = 0;
                        }
                }
        }
*/
//      fftwf_plan plan = fftwf_plan_r2r_1d(len, x, y, FFTW_RODFT11, FFTW_ESTIMATE);
        fftwf_plan plan = fftwf_plan_dft_r2c_3d(len, len, len, static_cast<float *>(input), y, FFTW_ESTIMATE);

        fftwf_execute(plan);

        for (int i = 0; i < len*len*lena; i++)
                z[i] = sqrt(y[i][REAL]*y[i][REAL]+y[i][IMAG]*y[i][IMAG]);

        for (int i = len-1; i > len/2; i-- ) {
                for (int j = 0; j < len; j++) {
                        for (int k = 1; k < lena; k++) {
                                z[(len-i)*len*lena+j*lena+k] = z[(len-i)*len*lena+j*lena+k] + z[i*len*lena+j*lena+k];
                        }
                }
        }

        for (int i = 0; i < len; i++ ) {
                for (int j = len-1; j > len/2; j--) {
                        for (int k = 1; k < lena; k++) {
                                z[i*len*lena+(len-j)*lena+k] = z[i*len*lena+(len-j)*lena+k] + z[i*len*lena+j*lena+k];
                        }
                }
        }

        fftwf_destroy_plan(plan);

	energy_low = 0;
	for (int i = 0; i < lena/2+1; i++)
		for (int j = 0; j < lena/2+1; j++)
			for (int k = 0; k < lena/2+1; k++)
				energy_low += z[i*len*lena+j*lena+k];
	energy_low = energy_low/pow(lena/2+1, 3);
	printf("%f ", energy_low);

	energy_mid_low = 0;
	for (int i = lena/2+1; i < lena; i++)
                for (int j = 0; j < lena/2+1; j++)
                        for (int k = 0; k < lena/2+1; k++)
				energy_mid_low += z[i*len*lena+j*lena+k];
        for (int i = 0; i < lena/2+1; i++)
                for (int j = 0; j < lena/2+1; j++)
                        for (int k = lena/2+1; k < lena; k++)
                                energy_mid_low += z[i*len*lena+j*lena+k];
        for (int i = 0; i < lena/2+1; i++)
                for (int j = lena/2+1; j < lena; j++)
                        for (int k = 0; k < lena/2+1; k++)
                                energy_mid_low += z[i*len*lena+j*lena+k];
	energy_mid_low = energy_mid_low/3/((lena/2)*(lena/2+1)*(lena/2+1));
	printf("%f ", energy_mid_low);

	energy_mid_high = 0;
	for (int i = 0; i < lena/2+1; i++)
                for (int j = lena/2+1; j < lena; j++)
                        for (int k = lena/2+1; k < lena; k++)
                                energy_mid_high += z[i*len*lena+j*lena+k];
	for (int i = lena/2+1; i < lena; i++)
                for (int j = lena/2+1; j < lena; j++)
                        for (int k = 0; k < lena/2+1; k++)
                                energy_mid_high += z[i*len*lena+j*lena+k];
        for (int i = lena/2+1; i < lena; i++)
                for (int j = 0; j < lena/2+1; j++)
                        for (int k = lena/2+1; k < lena; k++)
                                energy_mid_high += z[i*len*lena+j*lena+k];
	energy_mid_high = energy_mid_high/3/((lena/2)*(lena/2)*(lena/2+1));
	printf("%f ", energy_mid_high);

	energy_high = 0;
        for (int i = lena/2+1; i < lena; i++)
                for (int j = lena/2+1; j < lena; j++)
                        for (int k = lena/2+1; k < lena; k++)
                                energy_high += z[i*len*lena+j*lena+k];
	energy_high = energy_high/pow((lena/2), 3);
	printf("%f ", energy_high);

	fft3dClock.stop();
	double fft3d_time = fft3dClock.getDuration();
	printf("%f ", fft3d_time);


	printf("%f ", minmax);

        //free(pos);
//      abs = (1/entropy)/800;
//	abs = (double) ((double) abs * (double) entropy);
//	double absTol = 5.0;
	abs = (double) 5.0;

//	abs = strConvert::to_double( compressorParameters["abs"] ) / (double) 10;
	printf("%f ", abs);





    Timer cTime; 
    cTime.start();

	zfp_type type = getZfpType( dataType );

    // allocate meta data for the 1D input array
    zfp_field* field;
    if (numDims == 1)
        field = zfp_field_1d(input, type, numel);
    else if (numDims == 2)
        field = zfp_field_2d(input, type, n[1], n[0]);
    else if (numDims == 3)
        field = zfp_field_3d(input, type, n[2], n[1], n[0]);
    else if (numDims == 4)
        field = zfp_field_4d(input, type, n[3], n[2], n[1], n[0]);
   


    // allocate meta data for a compressed stream
    zfp_stream* zfp = zfp_stream_open(NULL);

    // set absolute/relative error tolerance
    if (compressionTypeAbs)
        zfp_stream_set_accuracy(zfp, abs);
    else
        zfp_stream_set_precision(zfp, rel);

   	//allocate buffer for compressed data
    size_t bufsize = zfp_stream_maximum_size(zfp, field);
    output = malloc(bufsize);

    // associate bit stream with allocated buffer
    bitstream* stream = stream_open(output, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);

    // Compress
    size_t zfpsize = zfp_compress(zfp, field);
    if (!zfpsize)
    {
        std::cout << "compression failed\n";
        return 0;
    }

    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);

    zfpCompressedSize = zfpsize;
    cbytes = zfpsize;

    cTime.stop();


    log << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << cbytes << ", cRatio: " << (dataTypeSize*numel / (float)cbytes) << ", #elements: " << numel << std::endl;
    log << compressorName << " ~ CompressTime: " << cTime.getDuration() << " s " << std::endl;

    return 1;
}


inline int ZFPCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
    int numDims = 1;
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
        {
			numel *= n[i];
            numDims++;
        }

    // Read in json compression parameters
    double abs = 1E-3;
    int rel = 32;
    bool compressionTypeAbs = true;
    if ( compressorParameters.find("abs") != compressorParameters.end() )
    {
        abs = strConvert::to_double( compressorParameters["abs"] );
        compressionTypeAbs = true;
    }
    else if ( compressorParameters.find("rel") != compressorParameters.end())
    {
        rel = strConvert::to_int( compressorParameters["rel"] );
        compressionTypeAbs = false;
    }

    abs = 5.0;


    Timer dTime; 
    dTime.start();
    
    zfp_type type = getZfpType( dataType );

    // allocate meta data  array of decompressed data
    output = malloc(numel*dataTypeSize);

    zfp_field* field;
    if (numDims == 1)
        field = zfp_field_1d(output, type, numel);
    else if (numDims == 2)
        field = zfp_field_2d(output, type, n[1], n[0]);
    else if (numDims == 3)
        field = zfp_field_3d(output, type, n[2], n[1], n[0]);
    else if (numDims == 4)
        field = zfp_field_4d(output, type, n[3], n[2], n[1], n[0]);



    // allocate meta data for a compressed stream
    zfp_stream* zfp = zfp_stream_open(NULL);
   
    // set absolute/relative error tolerance
    if (compressionTypeAbs)
        zfp_stream_set_accuracy(zfp, abs);
    else
        zfp_stream_set_precision(zfp, rel);

    // allocate buffer for decompressed data and transfer data
    size_t bufsize = zfp_stream_maximum_size(zfp, field);
    void  *buffer = malloc(bufsize);
    memcpy(buffer, input, zfpCompressedSize);

    // associate bit stream with allocated buffer
    bitstream* stream = stream_open(buffer, zfpCompressedSize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);


    // DeCompress
    size_t zfpsize  = zfp_decompress(zfp, field);
    if (! zfpsize)
    {
        std::cout << "decompression failed\n";
        return 0;
    }
    
    free(buffer);
    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);

    dTime.stop();

    log << compressorName << " ~ DecompressTime: " << dTime.getDuration() << " s " << std::endl;

    return 1;
}

inline void ZFPCompressor::close()
{

}

#endif // _ZFP_COMPRESSOR_H_
#endif // CBENCH_HAS_ZFP

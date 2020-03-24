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

#include <zfp.h>

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

    Timer clock("compress");

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

    clock.stop("compress");


    debugLog << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << cbytes << ", cRatio: " << (dataTypeSize*numel / (float)cbytes) << ", #elements: " << numel << std::endl;
    debugLog << compressorName << " ~ CompressTime: " << clock.getDuration("compress") << " s " << std::endl;

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


    Timer clock("decompress");
    
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

    clock.stop("decompress");

    debugLog << compressorName << " ~ DecompressTime: " << clock.getDuration("decompress") << " s " << std::endl;

    return 1;
}

inline void ZFPCompressor::close()
{

}

#endif // _ZFP_COMPRESSOR_H_
#endif // CBENCH_HAS_ZFP

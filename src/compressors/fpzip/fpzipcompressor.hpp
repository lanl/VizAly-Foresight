#ifdef CBENCH_HAS_FPZIP

#ifndef _FPZIP_COMPRESSOR_H_
#define _FPZIP_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
#include "timer.hpp"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "fpzip.h"

class FPZIPCompressor: public CompressorInterface
{
    
  public:
    FPZIPCompressor();
    ~FPZIPCompressor();

    void init();
    int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    void close();
};

inline FPZIPCompressor::FPZIPCompressor()
{
    compressorName = "fpzip";
}

inline FPZIPCompressor::~FPZIPCompressor()
{
    
}


inline void FPZIPCompressor::init()
{

}

inline int FPZIPCompressor::compress(void* input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	Timer cTime; 
    cTime.start();
  	//fprintf(stderr, "compressing to %s\n", "memory");
  	//double t = now();
    size_t numel = n[0];
    output = malloc(1024 + dataTypeSize*numel);

  	FPZ* fpz = fpzip_write_to_buffer(output, 1024 + dataTypeSize*numel);

  	if (dataType == "float")
  		fpz->type = FPZIP_TYPE_FLOAT;
  	else
  		fpz->type = FPZIP_TYPE_DOUBLE;

  	fpz->prec = 8;
  	fpz->nx = numel;
  	fpz->ny = 1;
  	fpz->nz = 1;
  	fpz->nf = 1;


  	// perform actual compression
  	size_t cbytes = fpzip_write(fpz, input);
  	if (!cbytes) 
  	{
    	fprintf(stderr, "compression failed: %s\n", fpzip_errstr[fpzip_errno]);
    	return 0;
  	}

  	fpzip_write_close(fpz);


  	//t = now() - t;
  	//fprintf(stderr, "in=%zu out=%zu ratio=%.2f seconds=%.3f MB/s=%.3f\n", inbytes, outbytes, (double)inbytes / outbytes, t, (double)inbytes / (1024 * 1024 * t));
  	cTime.stop();

    log << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << cbytes << ", cRatio: " << (dataTypeSize*numel / (float)cbytes) << std::endl;
    log << compressorName << " ~ CompressTime: " << cTime.getDuration() << " s " << std::endl;

    return 1;
}



inline int FPZIPCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
    Timer dTime; dTime.start();

	size_t numel = n[0];
    for (int i = 1; i < 5; i++)
		if (n[i] != 0)
            numel *= n[i];

	output = malloc(dataTypeSize*numel);

	FPZ* fpz = fpzip_read_from_buffer(input);
	if (dataType == "float")
  		fpz->type = FPZIP_TYPE_FLOAT;
  	else
  		fpz->type = FPZIP_TYPE_DOUBLE;

  	fpz->prec = 8;
  	fpz->nx = numel;
  	fpz->ny = 1;
  	fpz->nz = 1;
  	fpz->nf = 1;


	if (!fpzip_read(fpz, output)) 
	{
    	fprintf(stderr, "decompression failed: %s\n", fpzip_errstr[fpzip_errno]);
    	return EXIT_FAILURE;
  	}


  	fpzip_read_close(fpz);

	dTime.stop(); 

	std::free(input);	input=NULL;

	log << compressorName << " ~ DecompressTime: " << dTime.getDuration() << " s " << std::endl;

	return 1;
}

inline void FPZIPCompressor::close()
{
	//SZ_Finalize();
}

#endif // _SZ_COMPRESSOR_H_
#endif // CBENCH_HAS_SZ


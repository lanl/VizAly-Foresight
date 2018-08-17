#ifndef _BLOSC_COMPRESSOR_H_
#define _BLOSC_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
#include "timer.hpp"
#include "blosc.h"

class BLOSCCompressor: public CompressorInterface
{
    
  public:
    BLOSCCompressor();
    ~BLOSCCompressor();

    void init();
    int compress(void *input, void *&output, size_t dataType, size_t n);
    int decompress(void *&input, void *&output, size_t dataType, size_t n);
    void close();

	size_t cbytes;
};

inline BLOSCCompressor::BLOSCCompressor()
{
    compressorName = "BLOSC";
}

inline BLOSCCompressor::~BLOSCCompressor()
{
    
}


inline void BLOSCCompressor::init()
{
	blosc_init();
}

inline int BLOSCCompressor::compress(void *input, void *&output, size_t dataType, size_t n)
{
	// compress
	Timer cTime; cTime.start();
	// Default Input Params: {clevel=9, shuffle=1, sizeof(data), idatasize, input, output, odatasize);
	size_t isize = dataType*n;
	size_t osize = isize + BLOSC_MAX_OVERHEAD;

	output = std::malloc(isize); //byte array;
	osize = blosc_compress(9, 1, dataType, isize, input, output, osize);
	
	if (osize < 0)
	{
		throw std::runtime_error("Compression error. Error code: " + std::to_string(osize));
	}
	if (osize > 0)
	{
		output = std::realloc(output, osize);
	}
	cTime.stop();

	cbytes = osize;

	log << "\n" << compressorName << " ~ InputBytes: " << isize << ", OutputBytes: " << osize << ", cRatio: " << (isize/osize) << std::endl;
	log << compressorName << " ~ CompressTime: " << cTime.getDuration() << " s " << std::endl;

    return 1;
}

inline int BLOSCCompressor::decompress(void *&input, void *&output, size_t dataType, size_t n)
{
	Timer dTime; dTime.start();
	size_t osize = dataType*n;
	output = std::malloc(osize);
	size_t sz = blosc_decompress(input, output, osize);
	if (sz < 0)
		throw std::runtime_error("Decompression error. Error code: " + std::to_string(sz));
	
	std::free(input); input = NULL;

	dTime.stop(); log << compressorName << " ~ DecompressTime: " << dTime.getDuration() << " s " << std::endl;

    return 1;
}

inline void BLOSCCompressor::close()
{
	blosc_destroy();
}

#endif
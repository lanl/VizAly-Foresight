#ifdef CBENCH_HAS_BLOSC

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
    int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    void close();
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

inline int BLOSCCompressor::compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

	// compress
	Timer clock("compress");

	// Default Input Params: {clevel=9, shuffle=1, sizeof(data), idatasize, input, output, odatasize);

	size_t isize = dataTypeSize*numel;
	size_t osize = isize + BLOSC_MAX_OVERHEAD;

	output = std::malloc(isize); //byte array;
	osize = blosc_compress(9, 1, dataTypeSize, isize, input, output, osize);
	
	if (osize < 0)
	{
		throw std::runtime_error("Compression error. Error code: " + std::to_string(osize));
	}
	if (osize > 0)
	{
		output = std::realloc(output, osize);
	}
	clock.stop("compress");

	cbytes = osize;

	debugLog << "\n" << compressorName << " ~ InputBytes: " << isize << ", OutputBytes: " << osize << ", cRatio: " << (isize/(float)osize) << ", #elements: " << numel << std::endl;
	debugLog << compressorName << " ~ CompressTime: " << clock.getDuration("compress") << " s " << std::endl;

    return 1;
}

inline int BLOSCCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

	Timer clock("decompress");

	size_t osize = dataTypeSize*numel;
	output = std::malloc(osize);
	size_t sz = blosc_decompress(input, output, osize);
	if (sz < 0)
		throw std::runtime_error("Decompression error. Error code: " + std::to_string(sz));
	
	std::free(input); input = NULL;

	clock.stop("decompress"); 
	debugLog << compressorName << " ~ DecompressTime: " << clock.getDuration("decompress") << " s " << std::endl;

    return 1;
}

inline void BLOSCCompressor::close()
{
	blosc_destroy();
}

#endif //_BLOSC_COMPRESSOR_H_
#endif // CBENCH_HAS_BLOSC

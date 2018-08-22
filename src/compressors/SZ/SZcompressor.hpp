#ifndef _SZ_COMPRESSOR_H_
#define _SZ_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
#include "timer.hpp"
#include "sz.h"

class SZCompressor: public CompressorInterface
{
    
  public:
    SZCompressor();
    ~SZCompressor();

    void init();
    int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t n);
    int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t n);
    void close();

    size_t cbytes;
};

inline SZCompressor::SZCompressor()
{
    compressorName = "SZ";
}

inline SZCompressor::~SZCompressor()
{
    
}


inline void SZCompressor::init()
{

}

inline int SZCompressor::compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t n)
{
	Timer cTime; cTime.start();
	SZ_Init(NULL); 

	std::uint64_t csize = 0;
	std::uint8_t *cdata = SZ_compress_args(SZ_FLOAT, static_cast<float *>(input), &csize, PW_REL, 0, 0, 1E-3, 0, 0, 0, 0, n);

	output = cdata;
	cTime.stop();

	cbytes = csize;

	log << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*n << ", OutputBytes: " << csize << ", cRatio: " << (dataTypeSize*n / csize) << std::endl;
	log << compressorName << " ~ CompressTime: " << cTime.getDuration() << " s " << std::endl;

	return 1;
}

inline int SZCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t n)
{
	Timer dTime; dTime.start();
	output = SZ_decompress(SZ_FLOAT, static_cast<std::uint8_t *>(input), cbytes, 0, 0, 0, 0, n);
	dTime.stop(); 

	std::free(input);

	log << compressorName << " ~ DecompressTime: " << dTime.getDuration() << " s " << std::endl;

	return 1;
}

inline void SZCompressor::close()
{
	SZ_Finalize();
}

#endif

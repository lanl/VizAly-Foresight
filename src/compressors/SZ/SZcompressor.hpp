#ifdef CBENCH_HAS_SZ

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
    int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    void close();
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

inline int SZCompressor::compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

	Timer cTime; cTime.start();
	SZ_Init(NULL); 

	double relTol = 1E-3;
	std::unordered_map<std::string, std::string>::const_iterator got = compressorParameters.find("tolerance");
	if( got != compressorParameters.end() )
		if (compressorParameters["tolerance"] != "")
			relTol = strConvert::to_double( compressorParameters["tolerance"] );

        double lowAbs = 0.0;
        std::unordered_map<std::string, std::string>::const_iterator got = compressorParameters.find("lower");
        if( got != compressorParameters.end() )
                if (compressorParameters["lower"] != "")
                        lowAbs = strConvert::to_double( compressorParameters["lower"] );

        double upperAbs = 0.0;
        std::unordered_map<std::string, std::string>::const_iterator got = compressorParameters.find("upper");
        if( got != compressorParameters.end() )
                if (compressorParameters["upper"] != "")
                        upperAbs = strConvert::to_double( compressorParameters["upper"] );

	std::uint64_t csize = 0;
	std::uint8_t *cdata = SZ_compress_args(SZ_FLOAT, static_cast<float *>(input), &csize, PW_REL, lowAbs, upperAbs, relTol, n[4], n[3], n[2], n[1], n[0]);
	
	output = cdata;
	cTime.stop();

	cbytes = csize;

	log << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << csize << ", cRatio: " << (dataTypeSize*numel / (float)csize) << std::endl;
	log << compressorName << " ~ CompressTime: " << cTime.getDuration() << " s " << std::endl;

	return 1;
}

inline int SZCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

	Timer dTime; dTime.start();
	output = SZ_decompress(SZ_FLOAT, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	dTime.stop(); 

	std::free(input);	input=NULL;

	log << compressorName << " ~ DecompressTime: " << dTime.getDuration() << " s " << std::endl;

	return 1;
}

inline void SZCompressor::close()
{
	SZ_Finalize();
}

#endif // _SZ_COMPRESSOR_H_
#endif // CBENCH_HAS_SZ


#ifdef CBENCH_HAS_MGARD

#ifndef _MGARD_COMPRESSOR_H_
#define _MGARD_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
#include "timer.hpp"
#include "mgard_capi.h"

class MGARDCompressor: public CompressorInterface
{
    
  public:
    MGARDCompressor();
    ~MGARDCompressor();

    void init();
    int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    void close();
};

inline MGARDCompressor::MGARDCompressor()
{
    compressorName = "MGARD";
}

inline MGARDCompressor::~MGARDCompressor()
{
    
}


inline void MGARDCompressor::init()
{
	 
}

inline int MGARDCompressor::compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

	// Read in json compression parameters
	double tol = 1E-3;
	std::unordered_map<std::string, std::string>::const_iterator got = compressorParameters.find("tolerance");
	if( got != compressorParameters.end() )
		if (compressorParameters["tolerance"] != "")
			tol = strConvert::to_double( compressorParameters["tolerance"] );

	// Set compression parameters
	int iflag = 0; //0 -> float, 1 -> double

	if(dataTypeSize == 8)
		iflag = 1;

	int out_size;

	Timer cTime; cTime.start();

	// Create a copy of the input data because compressor will auto-destroy it
	void * in_buff = std::malloc(numel*dataTypeSize);
	memcpy (in_buff, input, numel*dataTypeSize);

	//mgard_compress(flag, in_buff, &out_size, nrow, ncol, nfib, &tol)
    output = mgard_compress(iflag, in_buff, &out_size, n[0], n[1], n[2], &tol );
	std::uint64_t csize = out_size;
	cbytes = csize;

	cTime.stop();

	log << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << csize << ", cRatio: " << (dataTypeSize*numel / (float)csize) << std::endl;
	log << compressorName << " ~ CompressTime: " << cTime.getDuration() << " s " << std::endl;

    return 1;
}

inline int MGARDCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

	// Read in json compression parameters
	double tol = 1E-3;
	std::unordered_map<std::string, std::string>::const_iterator got = compressorParameters.find("tolerance");
	if( got != compressorParameters.end() )
		if (compressorParameters["tolerance"] != "")
			tol = strConvert::to_double( compressorParameters["tolerance"] );

	// Set compression parameters
	int iflag = 0; //0 -> float, 1 -> double

	if(dataTypeSize == 8)
		iflag = 1;

	int out_size = cbytes;

	Timer dTime; dTime.start();

	output = mgard_decompress(iflag, static_cast<std::uint8_t *>(input), out_size, n[0], n[1], n[2] );

	std::uint64_t dsize = out_size; // Is out_size updated by mgard?

	dTime.stop();
	dTime.stop(); log << compressorName << " ~ DecompressTime: " << dTime.getDuration() << " s " << std::endl;

    return 1;
}

inline void MGARDCompressor::close()
{
	
}

#endif // _MGARD_COMPRESSOR_H_
#endif // CBENCH_HAS_MGARD

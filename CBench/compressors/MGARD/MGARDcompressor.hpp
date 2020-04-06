#ifdef CBENCH_HAS_MGARD

#ifndef _MGARD_COMPRESSOR_H_
#define _MGARD_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
#include "timer.hpp"
#include "mgard_api.h"

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

	if (dataTypeSize == 8)
		iflag = 1;

	int out_size[1];

	Timer clock("compress");

	// Create a copy of the input data because compressor will auto-destroy it
	void * in_buff = std::malloc(numel*dataTypeSize);
	memcpy(in_buff, input, numel*dataTypeSize);

	//mgard_compress(flag, in_buff, &out_size, nrow, ncol, nfib, &tol)
	// Note: "tol" must be of same "type" as set iflagz
	std::cout << "n[0]: " << n[0] << std::endl;
	std::cout << "n[0]: " << n[1] << std::endl;
	std::cout << "n[0]: " << n[2] << std::endl;

	int _n[3];
	_n[0] = n[0];
	_n[1] = n[1];
	_n[2] = n[2];

	float *_data = new float[numel];
	//_data = in_buff;
	unsigned char * compressed_data = mgard_compress(iflag, _data, out_size, _n[0], _n[1], _n[2], tol );
	std::uint64_t csize = out_size[0];
	cbytes = csize;

	clock.stop("compress");

	debugLog << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << csize << ", cRatio: " << (dataTypeSize*numel / (float)csize) << ", #elements: " << numel << std::endl;
	debugLog << compressorName << " ~ CompressTime: " << clock.getDuration("compress") << " s " << std::endl;

	return 1;
}

inline int MGARDCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

	// Set compression parameters
	int iflag = 0; //0 -> float, 1 -> double

	if(dataTypeSize == 8)
		iflag = 1;

	int out_size = cbytes;

	Timer clock("decompress");

	//output = mgard_decompress(iflag, static_cast<std::uint8_t *>(input), out_size, n[0], n[1], n[2] );

	std::uint64_t dsize = out_size; // Is out_size updated by mgard?

	clock.stop("decompress");
	debugLog << compressorName << " ~ DecompressTime: " << clock.getDuration("decompress") << " s " << std::endl;

	return 1;
}

inline void MGARDCompressor::close()
{
	
}

#endif // _MGARD_COMPRESSOR_H_
#endif // CBENCH_HAS_MGARD

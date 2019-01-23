#ifdef CBENCH_HAS_BIG_CRUNCH

#ifndef _BIGCRUNCH_COMPRESSOR_H_
#define _BIGCRUNCH_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
#include "timer.hpp"
#include "bigcrunch/bigcrunch.hpp"

class BigCrunchCompressor: public CompressorInterface
{
    
  public:
    BigCrunchCompressor();
    ~BigCrunchCompressor();

    void init();
    int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    void close();
};

inline BigCrunchCompressor::BigCrunchCompressor()
{
    compressorName = "BigCrunch";
}

inline BigCrunchCompressor::~BigCrunchCompressor()
{
    
}


inline void BigCrunchCompressor::init()
{
	 
}

inline int BigCrunchCompressor::compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
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

	// Convert tol for BigCrunch input
	int tolErr=0;
	if(tol < 1)
	{
		while(tol < 1)
		{
			tolErr -=1;
			tol*=10;
		}
	}

	// Default Params { Error:-3, Tolerance:1, BLOSC_NTHREADS: 1, BLOSCFILTER:SHUFFLE, BLOSC_COMPRESSOR:ZSTD }
	bigcrunch::setting_t settings = { {bigcrunch::config_t::ERR, tolErr},
				 {bigcrunch::config_t::TOLERANCE, 0},
				 {bigcrunch::config_t::BLOSC_NTHREADS, 1},
				 {bigcrunch::config_t::BLOSC_COMPRESSOR, bigcrunch::blosc_compressor_t::ZSTD},
				 {bigcrunch::config_t::BLOSC_CLEVEL, 9},
				 {bigcrunch::config_t::BLOSC_FILTER, bigcrunch::blosc_filter_t::SHUFFLE}
	};

	Timer cTime; cTime.start();

	bigcrunch::bigcrunch bc(settings);

	std::uint8_t *cdata = nullptr;
	std::uint64_t csize = bc.compress(bigcrunch::darray(static_cast<float *>(input), numel), &cdata);

	output = cdata;
	cTime.stop();

	cbytes = csize;

	log << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << csize << ", cRatio: " << (dataTypeSize*numel / (float)csize) << std::endl;
	log << compressorName << " ~ CompressTime: " << cTime.getDuration() << " s " << std::endl;

    return 1;
}

inline int BigCrunchCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	// Read in json compression parameters
	double tol = 1E-3;
	std::unordered_map<std::string, std::string>::const_iterator got = compressorParameters.find("tolerance");
	if( got != compressorParameters.end() )
		if (compressorParameters["tolerance"] != "")
			tol = strConvert::to_double( compressorParameters["tolerance"] );

	// Convert tol for BigCrunch input
	int tolErr=0;
	if(tol < 1)
	{
		while(tol < 1)
		{
			tolErr -=1;
			tol*=10;
		}
	}

	bigcrunch::setting_t settings = { {bigcrunch::config_t::ERR, tolErr},
				 {bigcrunch::config_t::TOLERANCE, 0},
				 {bigcrunch::config_t::BLOSC_NTHREADS, 1},
				 {bigcrunch::config_t::BLOSC_COMPRESSOR, bigcrunch::blosc_compressor_t::ZSTD},
				 {bigcrunch::config_t::BLOSC_CLEVEL, 9},
				 {bigcrunch::config_t::BLOSC_FILTER, bigcrunch::blosc_filter_t::SHUFFLE}
	};
	Timer dTime; dTime.start();
	bigcrunch::bigcrunch bc(settings);

	auto rdata_array = bc.decompress(static_cast<std::uint8_t *>(input), cbytes);

	std::free(input);  input = NULL;
	output = rdata_array.data<float>();

	dTime.stop();
	dTime.stop(); log << compressorName << " ~ DecompressTime: " << dTime.getDuration() << " s " << std::endl;

    return 1;
}

inline void BigCrunchCompressor::close()
{
	
}

#endif // _BIGCRUNCH_COMPRESSOR_H_
#endif // CBENCH_HAS_BIG_CRUNCH

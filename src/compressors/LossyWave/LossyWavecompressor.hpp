#ifdef CBENCH_HAS_LOSSY_WAVE

#ifndef _LOSSYWAVE_COMPRESSOR_H_
#define _LOSSYWAVE_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
#include "timer.hpp"
#include "lossywave.hpp"

class LossyWaveCompressor: public CompressorInterface
{
    
  public:
    LossyWaveCompressor();
    ~LossyWaveCompressor();

    void init();
    int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t n);
    int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t n);
    void close();
};

inline LossyWaveCompressor::LossyWaveCompressor()
{
    compressorName = "LossyWave";
}

inline LossyWaveCompressor::~LossyWaveCompressor()
{
    
}


inline void LossyWaveCompressor::init()
{
	 
}

inline int LossyWaveCompressor::compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t n)
{
	// Set compression parameters
	int args[13] = { 404, 0, 128, 0,
					n, n, n,
					n, n, n,
					dataTypeSize, 50, 0 };

	lossywave::lossywave lw(args);
	
	Timer cTime; cTime.start();

	std::uint64_t csize = lw.compress(input, dataTypeSize, output);

	cTime.stop();

	cbytes = csize+4; //4 byte header

	log << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*n << ", OutputBytes: " << csize << ", cRatio: " << (dataTypeSize*n / (float)csize) << std::endl;
	log << compressorName << " ~ CompressTime: " << cTime.getDuration() << " s " << std::endl;

    return 1;
}

inline int LossyWaveCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t n)
{
	// Set compression parameters
	int args[13] = { 404, 0, 128, 0,
					n, n, n,
					n, n, n,
					dataTypeSize, 50, 0 };

	lossywave::lossywave lw(args);

	Timer dTime; dTime.start();

	std::uint64_t dsize = lw.decompress(input, output);

	dTime.stop();
	dTime.stop(); log << compressorName << " ~ DecompressTime: " << dTime.getDuration() << " s " << std::endl;

    return 1;
}

inline void LossyWaveCompressor::close()
{
	
}

#endif // _LOSSYWAVE_COMPRESSOR_H_
#endif // CBENCH_HAS_LOSSY_WAVE

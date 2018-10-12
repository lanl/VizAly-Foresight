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
    int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
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

inline int LossyWaveCompressor::compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

	int quant = 0;
	quant = strConvert::to_int(compressorParameters["quant"]);
	int pcnt = 50;
	pcnt = strConvert::to_int(compressorParameters["pcnt"]);

	// Set compression parameters
	int args[13] = { 404, 0, 128+quant, 0,
					n[0], n[1], n[2],
					n[0], n[1], n[2],
					dataTypeSize, pcnt, 0 };

	lossywave::lossywave lw(args);
	
	Timer cTime; cTime.start();

	output = std::malloc(numel*dataTypeSize);
	std::uint64_t csize = lw.compress(input, dataTypeSize, output);

	cTime.stop();

	cbytes = csize+4; //add 4 byte header for lz4

	log << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << csize << ", cRatio: " << (dataTypeSize*numel / (float)csize) << std::endl;
	log << compressorName << " ~ CompressTime: " << cTime.getDuration() << " s " << std::endl;

    return 1;
}

inline int LossyWaveCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

	int quant = 0;
	quant = strConvert::to_int(compressorParameters["quant"]);
	int pcnt = 50;
	pcnt = strConvert::to_int(compressorParameters["pcnt"]);

	// Set compression parameters
	int args[13] = { 404, 0, 128+quant, 0,
					n[0], n[1], n[2],
					n[0], n[1], n[2],
					dataTypeSize, pcnt, 0 };

	lossywave::lossywave lw(args);

	Timer dTime; dTime.start();
	output = std::malloc(numel*dataTypeSize);
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

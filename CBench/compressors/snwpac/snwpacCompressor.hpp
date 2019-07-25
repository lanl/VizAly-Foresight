#ifdef CBENCH_HAS_SNWPAC

#ifndef _SNWPAC_COMPRESSOR_H_
#define _SNWPAC_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
#include "timer.hpp"
#include "snwpac.hpp"

class SnwPacCompressor : public CompressorInterface
{
    
  public:
	SnwPacCompressor();
    ~SnwPacCompressor();

    void init();
    int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    void close();
};

inline SnwPacCompressor::SnwPacCompressor()
{
    compressorName = "snwpac";
}

inline SnwPacCompressor::~SnwPacCompressor()
{
    
}


inline void SnwPacCompressor::init()
{
	 
}

inline int SnwPacCompressor::compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

    int quant = 0;
    std::unordered_map<std::string, std::string>::const_iterator got = compressorParameters.find("quant");
    if( got != compressorParameters.end() )
        if (compressorParameters["quant"] != "")
	        quant = strConvert::to_int(compressorParameters["quant"]);

	int pcnt = 30;
    got = compressorParameters.find("pcnt");
    if (got != compressorParameters.end())
        if (compressorParameters["pcnt"] != "")
            pcnt = strConvert::to_int(compressorParameters["pcnt"]);

	// Set compression parameters
	int args[13] = { 404, 0, 128+quant, 0,
                    (int)n[0], (int)n[1], (int)n[2],
                    (int)n[0], (int)n[1], (int)n[2],
                    (int)dataTypeSize, pcnt, 0 };

	snwpac::snwpac sp(args,false);

	Timer cTime; cTime.start();

	output = std::malloc(numel*dataTypeSize);
	std::uint64_t csize = sp.compress(input, dataTypeSize, output);

	cTime.stop();

	cbytes = csize+4; //add 4 byte header for lz4

	log << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << csize << ", cRatio: " << (dataTypeSize*numel / (float)csize) << ", #elements: " << numel << std::endl;
	log << compressorName << " ~ CompressTime: " << cTime.getDuration() << " s " << std::endl;

    return 1;
}

inline int SnwPacCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

	int quant = 0;
    std::unordered_map<std::string, std::string>::const_iterator got = compressorParameters.find("quant");
    if (got != compressorParameters.end())
        if (compressorParameters["quant"] != "")
            quant = strConvert::to_int(compressorParameters["quant"]);

	int pcnt = 30;
    got = compressorParameters.find("pcnt");
    if (got != compressorParameters.end())
        if (compressorParameters["pcnt"] != "")
            pcnt = strConvert::to_int(compressorParameters["pcnt"]);

	// Set compression parameters
	int args[13] = { 404, 0, 128+quant, 0,
                    (int)n[0], (int)n[1], (int)n[2],
                    (int)n[0], (int)n[1], (int)n[2],
                    (int)dataTypeSize, pcnt, 0 };

	snwpac::snwpac sp(args,false);

	Timer dTime; dTime.start();
	output = std::malloc(numel*dataTypeSize);
	std::uint64_t dsize = sp.decompress(input, output);

	dTime.stop();
	dTime.stop(); log << compressorName << " ~ DecompressTime: " << dTime.getDuration() << " s " << std::endl;

    return 1;
}

inline void SnwPacCompressor::close()
{
	
}

#endif // _SNWPAC_COMPRESSOR_H_
#endif // CBENCH_HAS_SNWPAC

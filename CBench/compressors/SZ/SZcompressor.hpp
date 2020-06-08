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

	int mode = PW_REL; // Default by Sheng, PW_REL = 10
	std::string _mode = "PW_REL";

	double relTol = 0.0;
	std::unordered_map<std::string, std::string>::const_iterator got = compressorParameters.find("rel");
	if( got != compressorParameters.end() )
		if (compressorParameters["rel"] != "")
		{
			relTol = strConvert::to_double(compressorParameters["rel"]);
			mode = REL;
			_mode = "REL";
		}

	double absTol = 0.0;
	got = compressorParameters.find("abs");
	if( got != compressorParameters.end() )
		if (compressorParameters["abs"] != "")
		{
			absTol = strConvert::to_double(compressorParameters["abs"]);
			mode = ABS;
			_mode = "ABS";
		}

	double powerTol = 0.0;
	got = compressorParameters.find("pw_rel");
	if( got != compressorParameters.end() )
		if (compressorParameters["pw_rel"] != "")
		{
			powerTol = strConvert::to_double(compressorParameters["pw_rel"]);
			mode = PW_REL;
			_mode = "PW_REL";
			// Unknown mode, just fill in input to SZ
		}


        float entropy = 0;
        float value_max = 0;
        float value_min = 3.402823466e+38F;
        int pos[128];
        for (int i = 0; i < 128; i++)
                pos[i] = 0;
        for (int i = 0; i < numel; i++) {
                if (static_cast<float *>(input)[i] > value_max)
                        value_max = static_cast<float *>(input)[i];
                if (static_cast<float *>(input)[i] < value_min)
                        value_min = static_cast<float *>(input)[i];
        }
        float minmax = value_max - value_min;
        printf("min = %f, max = %f, minmax = %f, numel = %d\n", value_min, value_max, minmax, numel);
        for (int i = 0; i < numel; i++)
                pos[(int)floor(((static_cast<float *>(input)[i] - value_min)/minmax) * 128)] += 1;
        for (int i = 0; i < 128; i++) {
                if (pos[i] != 0)
                        entropy += (-float(pos[i])/numel) * (std::log(float(pos[i])/numel)/std::log(2));
        }
        printf("entropy = %f ", entropy);
        //free(pos);
//      absTol = (1/entropy)/800;
        absTol = entropy * 10;
        printf("error bound = %f\n", absTol);


//      printf("data start = %f %f %f %f\n", static_cast<float *>(input)[0], static_cast<float *>(input)[1], static_cast<float *>(input)[134217726], static_cast<float *>(input)[134217727]);


	std::size_t csize = 0;
	std::uint8_t *cdata = SZ_compress_args(SZ_FLOAT, static_cast<float *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);

	output = cdata;
	cTime.stop();

	cbytes = csize;

	log << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << csize << ", cRatio: " << (dataTypeSize*numel / (float)csize) << ", #elements: " << numel << std::endl;
	log << " ~ Mode used: " << _mode << " abs: " << absTol << ", rel: " << relTol << ", pw_tol: " << powerTol << " val: " << n[4] << ", " << n[3] << ", " << n[2] << ", " <<n[1] << ", " << n[0] << std::endl;
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


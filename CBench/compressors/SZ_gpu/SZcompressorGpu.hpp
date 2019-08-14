#ifdef CBENCH_HAS_SZ_GPU

#ifndef _SZ_GPU_COMPRESSOR_H_
#define _SZ_GPU_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
#include "timer.hpp"
#include "sz.h"

class SZCompressorGpu: public CompressorInterface
{

  public:
    SZCompressorGpu();
    ~SZCompressorGpu();

    void init();
    int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    void close();
};

inline SZCompressorGpu::SZCompressorGpu()
{
    compressorName = "SZ_gpu";
}

inline SZCompressorGpu::~SZCompressorGpu()
{

}


inline void SZCompressorGpu::init()
{

}

inline int SZCompressorGpu::compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	sz_opencl_state* gpu_state;

    if(sz_opencl_init(&gpu_state) == SZ_NSCS) {
      throw std::runtime_error(std::string(sz_opencl_error_msg(gpu_state)));
    }

    auto gpu_compress = [&gpu_state](auto... args){
		return sz_compress_float3d_opencl(gpu_state, args...);
	};

	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

	Timer cTime; cTime.start();
	SZ_Init(NULL);

	int mode = ABS;
	std::string _mode = "ABS";

	double relTol = 0.0;
	std::unordered_map<std::string, std::string>::const_iterator got = compressorParameters.find("abs");
	// if( got != compressorParameters.end() )
	// 	if (compressorParameters["rel"] != "")
	// 	{
	// 		relTol = strConvert::to_double(compressorParameters["rel"]);
	// 		mode = REL;
	// 		_mode = "REL";
	// 	}

	double absTol = 0.0;
	//got = compressorParameters.find("abs");
	if( got != compressorParameters.end() )
		if (compressorParameters["abs"] != "")
		{
			absTol = strConvert::to_double(compressorParameters["abs"]);
			mode = ABS;
			_mode = "ABS";
		}

	double powerTol = 0.0;
	// got = compressorParameters.find("pw_rel");
	// if( got != compressorParameters.end() )
	// 	if (compressorParameters["pw_rel"] != "")
	// 	{
	// 		powerTol = strConvert::to_double(compressorParameters["pw_rel"]);
	// 		mode = PW_REL;
	// 		_mode = "PW_REL";
	// 		// Unknown mode, just fill in input to SZ
	// 	}

	float absEB = (float) absTol;

	std::size_t csize = 0;
	auto cdata = gpu_compress(static_cast<float *>(input), n[0], n[1], n[2], absEB, &csize);

	output = cdata;
	cTime.stop();

	cbytes = csize;

	log << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << csize << ", cRatio: " << (dataTypeSize*numel / (float)csize) << ", #elements: " << numel << std::endl;
	log << " ~ Mode used: " << _mode << " abs: " << absTol << ", rel: " << relTol << ", pw_tol: " << powerTol << " val: " << n[4] << ", " << n[3] << ", " << n[2] << ", " <<n[1] << ", " << n[0] << std::endl;
	log << compressorName << " ~ CompressTime: " << cTime.getDuration() << " s " << std::endl;

	return 1;
}


inline int SZCompressorGpu::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	sz_opencl_state* gpu_state;

    if(sz_opencl_init(&gpu_state) == SZ_NSCS) {
      throw std::runtime_error(std::string(sz_opencl_error_msg(gpu_state)));
    }

	auto gpu_decompress = [&gpu_state](auto... args){
		sz_decompress_float_opencl(gpu_state, args...);
	};

	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];

	Timer dTime; dTime.start();

	float* new_data;
	gpu_decompress(&new_data,
        /*r5*/0,/*r4*/0,/*r3*/n[2],/*r2*/n[1],/*r1*/n[0],
        /*s5*/0,/*s4*/0,/*s3*/0,         /*s2*/0,         /*s1*/0, /*start_positions*/
        /*e5*/0,/*e4*/0,/*s3*/n[2],/*e2*/n[1],/*e1*/n[0], /*end positions*/
        static_cast<std::uint8_t *>(input), cbytes);
	// output = SZ_decompress(SZ_FLOAT, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	dTime.stop();

	output = new_data;

	//std::free(input);

	input=NULL;

	log << compressorName << " ~ DecompressTime: " << dTime.getDuration() << " s " << std::endl;

	return 1;
}

inline void SZCompressorGpu::close()
{
	SZ_Finalize();
}

#endif // _SZ_GPU_COMPRESSOR_H_
#endif // CBENCH_HAS_SZ_GPU


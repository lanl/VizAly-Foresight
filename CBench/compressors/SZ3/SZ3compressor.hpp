#ifdef CBENCH_HAS_SZ3

#ifndef _SZ3_COMPRESSOR_H_
#define _SZ3_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
//#include "timer.hpp"
//#include "sz.h"
#include "SZ3/api/sz.hpp"

#define CONFIG_FILE "/projects/exasky/vis_compression/code/VizAly-Foresight/ExternalDependencies/SZ3/tools/sz3/sz3.config"

using namespace SZ;

class SZ3Compressor: public CompressorInterface
{

  public:
    SZ3Compressor();
    ~SZ3Compressor();

    void init();
    int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n);
    void close();
};

inline SZ3Compressor::SZ3Compressor()
{
    compressorName = "SZ3";
}

inline SZ3Compressor::~SZ3Compressor()
{

}


inline void SZ3Compressor::init()
{

}

inline int SZ3Compressor::compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	SZ::Config conf;

	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
			numel *= n[i];
	if (n[1] == 0)
		conf = SZ::Config(n[0]);
	else if (n[2] == 0)
		conf = SZ::Config(n[1], n[0]);
	else if (n[3] == 0)
		conf = SZ::Config(n[2], n[1], n[0]);
	else if (n[4] == 0)
		conf = SZ::Config(n[3], n[2], n[1], n[0]);
	else{
		conf = SZ::Config(n[4], n[3], n[2], n[1], n[0]);
	}

	conf.loadcfg(CONFIG_FILE);

	//Timer clock("compress");

	//int mode = PW_REL; // Default by Sheng, PW_REL = 10
	//std::string _mode = "PW_REL";

	double relTol = 0.0;
	std::unordered_map<std::string, std::string>::const_iterator got = compressorParameters.find("rel");
	if ( got != compressorParameters.end() )
		if (compressorParameters["rel"] != "")
		{
			relTol = strConvert::to_double(compressorParameters["rel"]);
			//mode = REL;
			//_mode = "REL";
			conf.errorBoundMode = SZ::EB_REL;
			conf.relErrorBound = relTol;
		}

	double absTol = 0.0;
	got = compressorParameters.find("abs");
	if ( got != compressorParameters.end() )
		if (compressorParameters["abs"] != "")
		{
			absTol = strConvert::to_double(compressorParameters["abs"]);
			//mode = ABS;
			//_mode = "ABS";
			conf.errorBoundMode = SZ::EB_ABS;
			conf.absErrorBound = absTol;
		}

	conf.cmprAlgo = SZ::ALGO_INTERP_LORENZO; //ALGO_INTERP_LORENZO(default) ALGO_INTERP, ALGO_LORENZO_REG


	size_t csize = 0;
	//std::uint8_t *cdata;
	char* cdata;
	if (dataType == "float")
		cdata = SZ_compress(conf, (float*) input, csize);
	 	//cdata = SZ_compress_args(SZ_FLOAT, static_cast<float *>(input),   &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	//else if (dataType == "double")
	//	cdata = SZ_compress_args(SZ_DOUBLE, static_cast<double *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	//else if (dataType == "Int64")
	//	cdata = SZ_compress_args(SZ_INT64, static_cast<int64_t *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	//else if (dataType == "Int32")
	//	cdata = SZ_compress_args(SZ_INT32, static_cast<int32_t *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	//else if (dataType == "Int16")
	//	cdata = SZ_compress_args(SZ_INT16, static_cast<int16_t *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	//else if (dataType == "Int8")
	//	cdata = SZ_compress_args(SZ_INT8, static_cast<int8_t *>(input),   &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	//else if (dataType == "UInt64")
	//	cdata = SZ_compress_args(SZ_UINT64, static_cast<uint64_t *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	//else if (dataType == "UInt32")
	//	cdata = SZ_compress_args(SZ_UINT32, static_cast<uint32_t *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	//else if (dataType == "UInt16")
	//	cdata = SZ_compress_args(SZ_UINT16, static_cast<uint16_t *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	//else if (dataType == "UInt8")
	//	cdata = SZ_compress_args(SZ_UINT8, static_cast<uint8_t *>(input),   &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	//else
	//{
	//	debugLog << "Data type: " << dataType << " not supported!" << std::endl;
	//	std::cout << "Data type: " << dataType << " not supported!" << std::endl;
	//	return -1;
	//}

	output = cdata;
	//clock.stop("compress");

	cbytes = csize;

	debugLog << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << csize << ", cRatio: " << (dataTypeSize*numel / (float)csize) << ", #elements: " << numel << std::endl;
	//debugLog << " ~ Mode used: " << _mode << " abs: " << absTol << ", rel: " << relTol << ", pw_tol: " << powerTol << " val: " << n[4] << ", " << n[3] << ", " << n[2] << ", " <<n[1] << ", " << n[0] << std::endl;
	//debugLog << compressorName << " ~ CompressTime: " << clock.getDuration("compress") << " s " << std::endl;

	return 1;
}


inline int SZ3Compressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	//size_t numel = n[0];
	//size_t for (int i = 1; i < 5; i++)
	//size_t 	if (n[i] != 0)
	//size_t 		numel *= n[i];

	//size_t Timer clock("decompress");

	//size_t if (dataType == "float")
	//size_t 	output = SZ_decompress(SZ_FLOAT,  static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	//size_t else if (dataType == "double")
	//size_t 	output = SZ_decompress(SZ_DOUBLE, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	//size_t else if (dataType == "Int64")
	//size_t 	output = SZ_decompress(SZ_INT64, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	//size_t else if (dataType == "Int32")
	//size_t 	output = SZ_decompress(SZ_INT32, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	//size_t else if (dataType == "Int16")
	//size_t 	output = SZ_decompress(SZ_INT16, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	//size_t else if (dataType == "Int8")
	//size_t 	output = SZ_decompress(SZ_INT8,  static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	//size_t else if (dataType == "UInt64")
	//size_t 	output = SZ_decompress(SZ_UINT64, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	//size_t else if (dataType == "UInt32")
	//size_t 	output = SZ_decompress(SZ_UINT32, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	//size_t else if (dataType == "UInt16")
	//size_t 	output = SZ_decompress(SZ_UINT16, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	//size_t else if (dataType == "UInt8")
	//size_t 	output = SZ_decompress(SZ_UINT8,  static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	//size_t else
	//size_t {
	//size_t 	debugLog << "Data type: " << dataType << " not supported!" << std::endl;
	//size_t 	std::cout << "Data type: " << dataType << " not supported!" << std::endl;
	//size_t 	return -1;
	//size_t }
	//size_t 

	//size_t clock.stop("decompress");

	//size_t std::free(input);	input=NULL;

	//size_t debugLog << compressorName << " ~ DecompressTime: " << clock.getDuration("decompress") << " s " << std::endl;

	return 1;
}

inline void SZ3Compressor::close()
{
	//SZ_Finalize();
}

#endif // _SZ_COMPRESSOR_H_
#endif // CBENCH_HAS_SZ

/*

abs - absErrBound refers to the absolute error bound, which is to limit the (de)compression
errors to be within an absolute error. For example, absErrBound=0.0001 means the
decompressed value must be in [V-0.0001,V+0.0001], where V is the original true
value. 

rel - relBoundRatio refers to value-range based relative bound ratio, which is to limit the
(de)compression errors by considering the global data value range size (i.e., taking into
account the range size (max_value - min_value)). For example, suppose
relBoundRatio is set to 0.01, and the data set is {100,101,102,103,104,...,110}. In this
case, the maximum value is 110 and the minimum is 100. So, the global value range
size is 110-100=10, and the error bound will be 10*0.01=0.1, from the perspective of
"relBoundRatio". 

pw_rel - pw_relBoundRatio refers to point-wise relative Bound Ratio. pw_relBountRatio is to
limit the (de)compression errors by considering the point-wise original data values. For
example, suppose pw_relBoundRatio is set to 0.01, and the data set is
{100,101,102,103,104,...,110}, so the compression errors will be limited to 
{1,1.01,1.02,....1.10} for the data points. This parameter is only valid when
errorBoundMode = PW_REL. 

*/




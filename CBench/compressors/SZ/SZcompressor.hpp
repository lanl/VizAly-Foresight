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

	Timer clock("compress");
	SZ_Init(NULL);

	int mode = PW_REL; // Default by Sheng, PW_REL = 10
	std::string _mode = "PW_REL";

	double relTol = 0.0;
	std::unordered_map<std::string, std::string>::const_iterator got = compressorParameters.find("rel");
	if ( got != compressorParameters.end() )
		if (compressorParameters["rel"] != "")
		{
			relTol = strConvert::to_double(compressorParameters["rel"]);
			mode = REL;
			_mode = "REL";
		}

	double absTol = 0.0;
	got = compressorParameters.find("abs");
	if ( got != compressorParameters.end() )
		if (compressorParameters["abs"] != "")
		{
			absTol = strConvert::to_double(compressorParameters["abs"]);
			mode = ABS;
			_mode = "ABS";
		}

	double powerTol = 0.0;
	got = compressorParameters.find("pw_rel");
	if ( got != compressorParameters.end() )
		if (compressorParameters["pw_rel"] != "")
		{
			powerTol = strConvert::to_double(compressorParameters["pw_rel"]);
			mode = PW_REL;
			_mode = "PW_REL";
		}

	std::size_t csize = 0;
	std::uint8_t *cdata;
	if (dataType == "float")
	 	cdata = SZ_compress_args(SZ_FLOAT, static_cast<float *>(input),   &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "double")
		cdata = SZ_compress_args(SZ_DOUBLE, static_cast<double *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "Int64")
		cdata = SZ_compress_args(SZ_INT64, static_cast<int64_t *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "Int32")
		cdata = SZ_compress_args(SZ_INT32, static_cast<int32_t *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "Int16")
		cdata = SZ_compress_args(SZ_INT16, static_cast<int16_t *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "Int8")
		cdata = SZ_compress_args(SZ_INT8, static_cast<int8_t *>(input),   &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "UInt64")
		cdata = SZ_compress_args(SZ_UINT64, static_cast<uint64_t *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "UInt32")
		cdata = SZ_compress_args(SZ_UINT32, static_cast<uint32_t *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "UInt16")
		cdata = SZ_compress_args(SZ_UINT16, static_cast<uint16_t *>(input), &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "UInt8")
		cdata = SZ_compress_args(SZ_UINT8, static_cast<uint8_t *>(input),   &csize, mode, absTol, relTol, powerTol, n[4], n[3], n[2], n[1], n[0]);
	else
	{
		debugLog << "Data type: " << dataType << " not supported!" << std::endl;
		std::cout << "Data type: " << dataType << " not supported!" << std::endl;
		return -1;
	}

	output = cdata;
	clock.stop("compress");

	cbytes = csize;

	debugLog << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << csize << ", cRatio: " << (dataTypeSize*numel / (float)csize) << ", #elements: " << numel << std::endl;
	debugLog << " ~ Mode used: " << _mode << " abs: " << absTol << ", rel: " << relTol << ", pw_tol: " << powerTol << " val: " << n[4] << ", " << n[3] << ", " << n[2] << ", " <<n[1] << ", " << n[0] << std::endl;
	debugLog << compressorName << " ~ CompressTime: " << clock.getDuration("compress") << " s " << std::endl;

	return 1;
}


inline int SZCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	Timer clock("decompress");

	if (dataType == "float")
		output = SZ_decompress(SZ_FLOAT,  static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "double")
		output = SZ_decompress(SZ_DOUBLE, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "Int64")
		output = SZ_decompress(SZ_INT64, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "Int32")
		output = SZ_decompress(SZ_INT32, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "Int16")
		output = SZ_decompress(SZ_INT16, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "Int8")
		output = SZ_decompress(SZ_INT8,  static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "UInt64")
		output = SZ_decompress(SZ_UINT64, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "UInt32")
		output = SZ_decompress(SZ_UINT32, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "UInt16")
		output = SZ_decompress(SZ_UINT16, static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	else if (dataType == "UInt8")
		output = SZ_decompress(SZ_UINT8,  static_cast<std::uint8_t *>(input), cbytes, n[4], n[3], n[2], n[1], n[0]);
	else
	{
		debugLog << "Data type: " << dataType << " not supported!" << std::endl;
		std::cout << "Data type: " << dataType << " not supported!" << std::endl;
		return -1;
	}
	

	clock.stop("decompress");

	std::free(input);	input=NULL;

	debugLog << compressorName << " ~ DecompressTime: " << clock.getDuration("decompress") << " s " << std::endl;

	return 1;
}

inline void SZCompressor::close()
{
	SZ_Finalize();
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




#ifdef CBENCH_HAS_MGARD

#ifndef _MGARD_COMPRESSOR_H_
#define _MGARD_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
#include "timer.hpp"
#include "compress_x.hpp"

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
	std::vector<mgard_x::SIZE> shape;
	shape.push_back(n[0]);

	int numDims = 1;
	size_t numel = n[0];
	for (int i = 1; i < 5; i++)
		if (n[i] != 0)
		{
			numel *= n[i];
			numDims++;
			shape.push_back(n[i]);
		}

	std::cout << "Shape" << shape[0] << ", " << shape[1] << ", " << shape[2] << std::endl;
	std::cout << "numel" << numel << std::endl;
	std::cout << "numDims" << numDims << std::endl;


	Timer clock("compress");


	// // Create a copy of the input data because compressor will auto-destroy it
	// void * in_buff = std::malloc(numel*dataTypeSize);
	// memcpy(in_buff, input, numel*dataTypeSize);
	
	mgard_x::data_type mDataType = mgard_x::data_type::Float;
	if (dataType == "double")
	{
		mDataType = mgard_x::data_type::Double;
		std::cout << "double" << std::endl;
	}
	else
		if (dataType == "float")
		{
			mDataType = mgard_x::data_type::Float;
			std::cout << "Float" << std::endl;
		}
		else
		{
			std::cout << "double not supported!!!" << std::endl;
			return 0;
		}

	mgard_x::error_bound_type errorType = mgard_x::error_bound_type::REL;
	double tol = 1E-3;
	if ( compressorParameters.find("abs") != compressorParameters.end() )
	{
		errorType = mgard_x::error_bound_type::ABS;
		tol = strConvert::to_double(compressorParameters["abs"]);
		std::cout << "Abs " << " : " << tol << std::endl;
	}
	else if ( compressorParameters.find("rel") != compressorParameters.end() )
	{
		errorType = mgard_x::error_bound_type::REL;
		tol = strConvert::to_double(compressorParameters["rel"]);
		std::cout << "Rel " << " : " << tol << std::endl;
	}
	else
	{
		std::cout << "error type not supported!!!" << std::endl;
		return 0;
	}


	mgard_x::Config config;
	config.lossless = mgard_x::lossless_type::Huffman_Zstd;
	config.dev_type = mgard_x::device_type::Auto;


	mgard_x::compress(numDims, 
					mDataType, 
					shape, 
					tol, 
					0.0,	// smoothness parameter
                    errorType, 
					input,
                    output, 
					cbytes, 
					config, 
					false);


	clock.stop("compress");

	debugLog << "\n" << compressorName << " ~ InputBytes: " << dataTypeSize*numel << ", OutputBytes: " << cbytes << ", cRatio: " << (dataTypeSize*numel / (float)cbytes) << ", #elements: " << numel << std::endl;
	debugLog << compressorName << " ~ CompressTime: " << clock.getDuration("compress") << " s " << std::endl;

	return 1;
}


inline int MGARDCompressor::decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n)
{
	Timer clock("decompress");

	mgard_x::Config config;
	config.lossless = mgard_x::lossless_type::Huffman_Zstd;
	config.dev_type = mgard_x::device_type::Auto;

	mgard_x::decompress(input, 
						cbytes,
                        output, 
						config, 
						false);

	clock.stop("decompress");
	debugLog << compressorName << " ~ DecompressTime: " << clock.getDuration("decompress") << " s " << std::endl;

	std::free(input);	input=NULL;


	return 1;
}

inline void MGARDCompressor::close()
{
	
}

#endif // _MGARD_COMPRESSOR_H_
#endif // CBENCH_HAS_MGARD

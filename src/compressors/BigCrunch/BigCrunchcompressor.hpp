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
    int compress(void *input, void *&output, size_t dataType, size_t n);
    int decompress(void *&input, void *&output, size_t dataType, size_t n);
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

inline int BigCrunchCompressor::compress(void *input, void *&output, size_t dataType, size_t n)
{
	
    return 1;
}

inline int BigCrunchCompressor::decompress(void *&input, void *&output, size_t dataType, size_t n)
{
	
    return 1;
}

inline void BigCrunchCompressor::close()
{
	
}

#endif
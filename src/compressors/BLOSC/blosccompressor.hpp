#ifndef _BLOSC_COMPRESSOR_H_
#define _BLOSC_COMPRESSOR_H_

#include <sstream>
#include "compressorInterface.hpp"
#include "timer.hpp"

class BLOSCCompressor: public CompressorInterface
{
    
  public:
    BLOSCCompressor();
    ~BLOSCCompressor();

    void init();
    int compress(void *input, void *output, size_t dataType, size_t n);
    int decompress(void *input, void *output, size_t dataType, size_t n);
    int close();
};

inline BLOSCCompressor::BLOSCCompressor()
{
    compressorName = "BLOSC";
}

inline BLOSCCompressor::~BLOSCCompressor()
{
    
}


inline void BLOSCCompressor::init()
{

}

inline int BLOSCCompressor::compress(void *input, void *output, size_t dataType, size_t n)
{

    return 1;
}

inline int BLOSCCompressor::decompress(void *input, void *output, size_t dataType, size_t n)
{

    return 1;
}

inline void BLOSCCompressor::close()
{


}

#endif
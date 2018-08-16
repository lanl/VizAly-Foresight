#ifndef _COMPRESSOR_INTERFACE_H_
#define _COMPRESSOR_INTERFACE_H_

#include <string>
#include <sstream>

class CompressorInterface
{
  protected:
    std::string compressorName;
    std::stringstream log;

  public:
    virtual void init() = 0;
    virtual int compress(void *input, void *output, size_t dataType, size_t n) = 0;
    virtual int decompress(void *input, void *output, size_t dataType, size_t n) = 0;
    virtual void close() = 0;

    std::string getCompressorInfo();
    std::string getLog() { return log.str(); }
};


inline std::string CompressorInterface::getCompressorInfo()
{
    std::stringstream dataInfo;
    dataInfo << "\nCompressor: " << compressorName << std::endl;

    return dataInfo.str();
}

#endif

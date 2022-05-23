/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/

#ifndef _COMPRESSOR_INTERFACE_H_
#define _COMPRESSOR_INTERFACE_H_

#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>

#include "log.hpp"
#include "timer.hpp"
#include "memory.hpp"
#include "strConvert.hpp"


class CompressorInterface
{
  public:
    std::unordered_map<std::string, std::string> compressorParameters;

  protected:
    std::string compressorName;		// internal compressor name
    size_t cbytes;					// compressed stream size in bytes

  public:
    virtual void init() = 0;
    virtual int compress(void *input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n) = 0;
    virtual int decompress(void *&input, void *&output, std::string dataType, size_t dataTypeSize, size_t * n) = 0;
    virtual void close() = 0;

    std::string getCompressorInfo();
    std::string getCompressorName(){ return compressorName; }
    size_t getCompressedSize(){ return cbytes; }
    void setCompressedSize(size_t _cbytes){ cbytes = _cbytes; }
    std::string getParamsInfo();
};


inline std::string CompressorInterface::getCompressorInfo()
{
    std::stringstream dataInfo;
    dataInfo << "\nCompressor: " << compressorName << std::endl;

    return dataInfo.str();
}


std::string CompressorInterface::getParamsInfo()
    {
        std::string paramString = "";
        for (auto it=compressorParameters.begin(); it!=compressorParameters.end(); it++)
        {
            if (paramString != "")
                paramString += "_";
            paramString += (*it).first + (*it).second;
        }

        return paramString;  
    }
#endif

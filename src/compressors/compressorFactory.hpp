/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/

#ifndef _COMPRESSOR_FACTORY_H_
#define _COMPRESSOR_FACTORY_H_

#include "compressorIncludes.h"
#include "compressorInterface.hpp"

class CompressorFactory
{
  public:
	static CompressorInterface * createCompressor(std::string compressorName)
	{
		if (compressorName == "blosc")
		  return new BLOSCCompressor();
	  #ifdef CBENCH_HAS_BIG_CRUNCH
		else if (compressorName == "BigCrunch")
		  return new BigCrunchCompressor();
	  #endif
	  #ifdef CBENCH_HAS_ZFP
		else if (compressorName == "zfp")
		 return new ZFPCompressor();
	  #endif
	  #ifdef CBENCH_HAS_SZ
		else if (compressorName == "SZ")
		  return new SZCompressor();
	  #endif
	  #ifdef CBENCH_HAS_LOSSY_WAVE
		else if (compressorName == "LossyWave")
		  return new LossyWaveCompressor();
	  #endif
	  #ifdef CBENCH_HAS_MGARD
		else if (compressorName == "MGARD")
		  return new MGARDCompressor();
	  #endif
	  #ifdef CBENCH_HAS_ISABELA
		else if (compressorName == "Isabela")
		  return new IsabelaCompressor();
	  #endif
	  #ifdef CBENCH_HAS_FPZIP
		else if (compressorName == "fpzip")
			return new FPZIPCompressor();
	  #endif
		else
		  return NULL;
	}
};


#endif

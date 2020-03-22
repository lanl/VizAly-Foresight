/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/

#ifndef _LOADER_FACTORY_H_
#define _LOADER_FACTORY_H_

#include "dataLoaderIncludes.h"
#include "dataLoaderInterface.hpp"

class DataLoaderFactory
{
  public:
	static DataLoaderInterface * createLoader(std::string loaderName)
	{
		if (loaderName == "HACC")
			return new HACCDataLoader();

	  #ifdef CBENCH_HAS_NYX
		if (loaderName == "NYX")
			return new NYXDataLoader();
	  #endif

	  #ifdef CBENCH_HAS_VTK
		if (loaderName == "VTI")
			return new VTKDataLoader();
	  #endif

	  #ifdef CBENCH_HAS_GDA
		if (loaderName == "GDA")
			return new GDADataLoader();
	  #endif

	  #ifdef CBENCH_HAS_RAW
		if (loaderName == "RAW")
			return new GenericBinaryLoader();
	  #endif

		return NULL;
	}
};

#endif
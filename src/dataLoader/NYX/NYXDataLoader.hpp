/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Jesus Pulido
================================================================================*/

#ifndef _NYX_LOADER_H_
#define _NYX_LOADER_H_

#include <sstream>
#include "dataLoaderInterface.hpp"
#include "timer.hpp"

class NYXDataLoader: public DataLoaderInterface
{
	int numRanks;
	int myRank;

  public:
	NYXDataLoader();
	~NYXDataLoader();
	int allocateMem(std::string dataType, size_t numElements, int offset);
	int deAllocateMem();

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int close() { return deAllocateMem(); }
};


inline NYXDataLoader::NYXDataLoader()
{
	myRank = 0;
	numRanks = 0;
	loader = "NYX";
}

inline NYXDataLoader::~NYXDataLoader()
{
	deAllocateMem();
}


inline void NYXDataLoader::init(std::string _filename, MPI_Comm _comm)
{
	filename = _filename;
	comm = _comm;

	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);
}



inline int NYXDataLoader::allocateMem(std::string dataType, size_t numElements, int offset)
{
	
	return 1;
}


inline int NYXDataLoader::deAllocateMem()
{
	if (data == NULL) // already deallocated!
		return 1;

	data = NULL;

	return 1;
}


inline int NYXDataLoader::loadData(std::string paramName)
{
	Timer clock;
	log.str("");

	clock.start();
	
	clock.stop();
	log << "Loading data took " << clock.getDuration() << " s" << std::endl;

	return 1; // All good
}

#endif
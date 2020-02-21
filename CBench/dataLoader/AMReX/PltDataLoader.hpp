/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/

#pragma once

#include <sstream>
#include <string>
#include <unordered_map>

#include "dataLoaderInterface.hpp"

#include "strConvert.hpp"
#include "json.hpp"
#include "timer.hpp"
#include "utils.hpp"

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_DataServices.H"  
#include "AMReX_Utility.H"
#include "AMReX_FArrayBox.H"

class PltDataLoader: public DataLoaderInterface
{
	int numRanks;
	int myRank;

	float origRealDims[3];

	void *tempData;

	amrex::ParmParse pp;

  public:
	PltDataLoader();
	~PltDataLoader();

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int saveCompData(std::string paramName, void * cData);
	int writeData(std::string _filename);

	int saveInputFileParameters() { return 1; };
	int close() { return deAllocateMem(dataType, data); }
	void setParam(std::string paramName, std::string type, std::string value) {};
	bool loadUncompressedFields(nlohmann::json const&) { return false; }
};


inline PltDataLoader::PltDataLoader()
{

}


inline PltDataLoader::~PltDataLoader()
{

}


inline void PltDataLoader::init(std::string _filename, MPI_Comm _comm)
{
	filename = _filename;
	comm = _comm;
	saveData = false;

	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);

	pp.get("input_path", filename);
	

	log << "original dims: " << origDims[0] << ", " << origDims[1] << ", " << origDims[1] << std::endl;
	log << "real dims : " << origRealDims[0] << ", " << origRealDims[1] << ", " << origRealDims[1] << std::endl;
	log << "dataType : " << dataType << std::endl;
}



inline int PltDataLoader::loadData(std::string paramName)
{

}



inline int PltDataLoader::saveCompData(std::string paramName, void * cData)
{
	Timer clock;
	clock.start();
	
	allocateMem(dataType, numElements, 0, tempData);
	memcpy(tempData, cData, numElements*getDataypeSize(dataType));

	clock.stop();
	log << "saving data took: " << clock.getDuration() << " s" << std::endl;
}



inline int PltDataLoader::writeData(std::string _filename)
{
	Timer clock;
	clock.start();

	
	clock.stop();
	log << "writing data took: " << clock.getDuration() << " s" << std::endl;
}



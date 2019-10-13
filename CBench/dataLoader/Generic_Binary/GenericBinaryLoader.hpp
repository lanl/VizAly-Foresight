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

struct scalar
{
	std::string name;
	size_t offset;
	std::string type;
	void *data;

	scalar(){};
	scalar(std::string _name, size_t _offset, std::string _type)
	{
		name = _name;
		offset = _offset;
		type = _type;
	}

	std::string toString()
	{
		std::string temp;
		temp = strConvert::toStr(name) + " " + strConvert::toStr(offset) + " " + strConvert::toStr(type);
		return temp;
	}
};


class GenericBinaryLoader: public DataLoaderInterface
{
	int numRanks;
	int myRank;

	float origRealDims[3];
	MPI_Datatype mpiDataType;
	int mpiDivisions[3];

	std::vector<scalar> scalars;

  public:
	GenericBinaryLoader();
	~GenericBinaryLoader();

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int saveCompData(std::string paramName, void * cData);
	int writeData(std::string _filename);

	int saveInputFileParameters() { return 1; };
	int close() { return deAllocateMem(dataType, data); }
	void setParam(std::string paramName, std::string type, std::string value) {};
	bool loadUncompressedFields(nlohmann::json const&) { return false; }
};


inline GenericBinaryLoader::GenericBinaryLoader()
{
	timestep = -1;
}


inline GenericBinaryLoader::~GenericBinaryLoader()
{

}


inline void GenericBinaryLoader::init(std::string _filename, MPI_Comm _comm)
{
	filename = _filename;
	comm = _comm;
	saveData = false;

	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);

	if (timestep != -1)
		filename = filename + "_" + strConvert::toStr(timestep);
	std::string metadataFile = filename + ".info";

	std::ifstream ifs(metadataFile.c_str());
	if (ifs.is_open()) 
	{
		ifs >> origDims[0] >> origDims[1] >> origDims[2];

		while (ifs.good()) 
   	 	{
   	 		scalar temp;
   	 		ifs >> temp.name >> temp.offset >> temp.type;
   	 		scalars.push_back(temp);
    	}
  	}
  	else 
  	{
    	if (myRank == 0)
			std::cout << "Unable to open file " << metadataFile << ". This program will now exit!" << std::endl;

		MPI_Finalize();
  	}



	log << "original dims: " << origDims[0] << ", " << origDims[1] << ", " << origDims[1] << std::endl;
	for (auto s: scalars)
  		log << s.toString() << std::endl;

	log << "dataType : " << dataType << std::endl;
}



inline int GenericBinaryLoader::loadData(std::string paramName)
{
	Timer clock;
	clock.start();

	float minVal, maxVal, avg;
	totalNumberOfElements = origDims[0] * origDims[1] * origDims[2];

	
	// MPI split!!!
	Partition current = getPartition(myRank, numRanks, origDims[0], origDims[1], origDims[2]);
	rankOffset[0] = current.min_x;
	rankOffset[1] = current.min_y;
	rankOffset[2] = current.min_z;

	sizePerDim[0] = current.max_x - current.min_x;
	sizePerDim[1] = current.max_y - current.min_y;
	sizePerDim[2] = current.max_z - current.min_z;

	numElements = sizePerDim[0] * sizePerDim[1] * sizePerDim[2];


	// Get scalar parameters
	std::string currentDataType;
	size_t offset;
	bool found = false;
	int index = 0;
	for (auto s: scalars)
		if (paramName == s.name)
		{
			found = true;
			currentDataType = s.type;
			offset = s.offset;
			index++;
		}

	if (found == false)
		return -1;


	// Read the binary file
	if (timestep != -1)
		filename = filename + "_" + strConvert::toStr(timestep);
	std::string dataFile = filename + ".raw";


	// MPI collective File Read
	MPI_File fh;
	MPI_Status status;
	MPI_Datatype filetype;
	mpiDataType = getMPIType(dataType);

	int distribs[3];
	distribs[0] = MPI_DISTRIBUTE_BLOCK;;
	distribs[1] = MPI_DISTRIBUTE_BLOCK;;
	distribs[2] = MPI_DISTRIBUTE_BLOCK;;

	int dargs[3];
	dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
	dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
	dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;

	getMPIDivisions(numRanks, 3, mpiDivisions);


	int fullsize[3];
	for (int i = 0; i < 3; i++)
		fullsize[i] = origDims[i];

	MPI_Type_create_darray(numRanks, myRank, 3, fullsize, distribs, dargs, mpiDivisions, MPI_ORDER_C, mpiDataType, &filetype);
	MPI_Type_commit(&filetype);

	MPI_File_open(comm, dataFile.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

	MPI_File_set_view(fh, 0, mpiDataType, filetype, "native", MPI_INFO_NULL);


	// Create space for data and store data there
	allocateMem(dataType, numElements, 0, scalars[index].data);

	if (dataType == "float")
       	MPI_File_read_at_all(fh, offset, ((float *)data), numElements, mpiDataType, &status);
    else if (dataType == "double")
        MPI_File_read_at_all(fh, offset, ((float *)data), numElements, mpiDataType, &status);
    else if (dataType == "int")
        MPI_File_read_at_all(fh, offset, ((float *)data), numElements, mpiDataType, &status);
    else
    {	
    	if (myRank == 0)
    		std::cout << "Unsupported data type " << dataType << " . The program will now exit!" << std::endl;

    	MPI_Finalize();
    }

	MPI_File_close( &fh );

	

	clock.stop();
	log << "Param: " << paramName << std::endl;
	log << "min: " << minVal << ", max: " << maxVal<< ", avg: " << avg << std::endl;
	log << "mpiDivisions: " << mpiDivisions[0] << ", " << mpiDivisions[1] << ", " << mpiDivisions[2] << std::endl;
	log << "sizePerDim: "	<< sizePerDim[0]	<< ", " << sizePerDim[1]<< ", " << sizePerDim[2]<< std::endl;
	log << "rankOffset: "	<< rankOffset[0]	<< ", " << rankOffset[1]<< ", " << rankOffset[2]<< std::endl;
	log << "numElements: "	<< numElements << std::endl;
	log << "totalNumberOfElements: " << totalNumberOfElements << std::endl;

	log << "Loading data took: " << clock.getDuration() << " s" << std::endl;
}



inline int GenericBinaryLoader::saveCompData(std::string paramName, void * cData)
{
	Timer clock;
	clock.start();
	
	for (auto s: scalars)
		if (paramName == s.name)
			memcpy(s.data, cData, numElements*getDataypeSize(s.type));

	clock.stop();
	log << "saving data took: " << clock.getDuration() << " s" << std::endl;
}



inline int GenericBinaryLoader::writeData(std::string _filename)
{
	Timer clock;
	clock.start();

	MPI_File fh;
	MPI_Status status;
	MPI_Datatype filetype;

	int fullsize[3];
	for (int i = 0; i < 3; i++)
		fullsize[i] = origDims[i];

	int distribs[3];
	distribs[0] = MPI_DISTRIBUTE_BLOCK;;
	distribs[1] = MPI_DISTRIBUTE_BLOCK;;
	distribs[2] = MPI_DISTRIBUTE_BLOCK;;

	int dargs[3];
	dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
	dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
	dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;

	MPI_Type_create_darray(numRanks, myRank, 3, fullsize, distribs, dargs, mpiDivisions, MPI_ORDER_C, mpiDataType, &filetype);
	MPI_Type_commit(&filetype);

	_filename = _filename + ".raw";
	MPI_File_open(comm, _filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	MPI_File_set_view(fh, 0, mpiDataType, filetype, "native", MPI_INFO_NULL);
	
	for (auto s: scalars)
	{
		if (dataType == "float")
	       	MPI_File_write_at_all(fh, s.offset, (float *)s.data, numElements, mpiDataType, &status);
	    else if (dataType == "double")
	        MPI_File_write_at_all(fh, s.offset, (double *)s.data, numElements, mpiDataType, &status);
	    else if (dataType == "int")
	        MPI_File_write_at_all(fh, s.offset, (int *)s.data, numElements, mpiDataType, &status);

	    deAllocateMem(s.type, s.data);
	}

	MPI_File_close( &fh );

	
	clock.stop();
	log << "writing data took: " << clock.getDuration() << " s" << std::endl;
}

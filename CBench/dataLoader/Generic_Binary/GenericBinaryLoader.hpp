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
	bool toBeCompressed;
	bool compressed;
	void *data;

	scalar()
	{
		toBeCompressed = true;
		compressed = false;
		data = NULL;
	}

	scalar(std::string _name, size_t _offset, std::string _type)
	{
		name = _name;
		offset = _offset;
		type = _type;
		size_t numLocalElements;
		toBeCompressed = true;
		compressed = false;
		data = NULL;
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

	std::string rawDataFileName;

	std::vector<scalar> scalars;

  public:
	GenericBinaryLoader();
	~GenericBinaryLoader();

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int saveCompData(std::string paramName, void * cData);
	int writeData(std::string _filename);

	int saveInputFileParameters(){ return 1; }
	int close() { return deAllocateMem(dataType, data); }
	void setParam(std::string paramName, std::string type, std::string value) {};
	bool loadUncompressedFields(nlohmann::json const&) { return false; }
};


inline GenericBinaryLoader::GenericBinaryLoader()
{
	timestep = -1;
	myRank = 0;
	numRanks = 0;
	loader = "RAW";
	saveData = false;
	data = NULL;
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
;

	if (timestep != -1)
		filename = filename + "_" + strConvert::toStr(timestep);
	std::string metadataFile = filename + ".info";

	std::ifstream ifs(metadataFile.c_str());
	if (ifs.is_open()) 
	{
		ifs >> rawDataFileName;
		ifs >> origDims[0] >> origDims[1] >> origDims[2];

		while (ifs.good()) 
   	 	{
   	 		scalar temp;
   	 		ifs >> temp.name >> temp.type >> temp.offset;
			if (temp.name != "")
   	 			scalars.push_back(temp);
    	}
  	}
  	else 
  	{
    	if (myRank == 0)
			std::cout << "Unable to open file " << metadataFile << ". This program will now exit!" << std::endl;

		MPI_Finalize();
  	}



  	log << "\nRaw data file: " << rawDataFileName << std::endl;
	log << "original dims: " << origDims[0] << ", " << origDims[1] << ", " << origDims[1] << std::endl;
	for (auto s: scalars)
  		log << s.toString() << std::endl;
  	log << std::endl;
}



inline int GenericBinaryLoader::loadData(std::string paramName)
{
	Timer clock;
	clock.start();


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

	int index;
	for (int i=0; i<scalars.size(); i++)
		if (paramName == scalars[i].name)
		{
			found = true;
			currentDataType = scalars[i].type;
			offset = scalars[i].offset;

			index = i;
			break;
		}

	if (found == false)
		return -1;


	elemSize = Memory::sizeOf[currentDataType];
	dataType = currentDataType;
	MPI_Offset mpi_offset = offset;


	//
	// Read the binary file

	// MPI collective File Read
	MPI_File fh;
	MPI_Status status;
	MPI_Datatype filetype;
	mpiDataType = getMPIType(currentDataType);
	std::stringstream minMaxAvgLog;

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

	MPI_Barrier(comm);

	MPI_Type_create_darray(numRanks, myRank, 3, fullsize, distribs, dargs, mpiDivisions, MPI_ORDER_C, mpiDataType, &filetype);
	MPI_Type_commit(&filetype);

	int rc = MPI_File_open(comm, rawDataFileName.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	if (rc) 
	{
		if (myRank == 0)
        	std::cout << "Unable to open file " << rawDataFileName << std::endl;
        return -1;
    }

	MPI_File_set_view(fh, mpi_offset, MPI_BYTE, filetype, "native", MPI_INFO_NULL);


	// Create space for data and store data there
	if (scalars[index].toBeCompressed)
	{
		allocateMem(currentDataType, numElements, 0, data);

		if (currentDataType == "float")
		{
	       	MPI_File_read_all(fh, ((float *)data), numElements, mpiDataType, &status);

	       	float minVal, maxVal, avg;
	        minMax( ((float *)data), numElements, minVal, maxVal, avg);
	        minMaxAvgLog << "min: " << minVal << ", max: " << maxVal<< ", avg: " << avg << std::endl;
		}
	    else if (currentDataType == "double")
	    {
	        MPI_File_read_all(fh, ((double *)data), numElements, mpiDataType, &status);

	        double minVal, maxVal, avg;
	        minMax( ((double *)data), numElements, minVal, maxVal, avg);
	        minMaxAvgLog << "min: " << minVal << ", max: " << maxVal<< ", avg: " << avg << std::endl;
	    }
	    else if (currentDataType == "int")
	    {
	        MPI_File_read_all(fh, ((int *)data), numElements, mpiDataType, &status);

	        int minVal, maxVal, avg;
	        minMax( ((int *)data), numElements, minVal, maxVal, avg);
	        minMaxAvgLog << "min: " << minVal << ", max: " << maxVal<< ", avg: " << avg << std::endl;
	    }
	    else
	    {	
	    	if (myRank == 0)
	    		std::cout << "Unsupported data type " << currentDataType << " . The program will now exit!" << std::endl;

	    	MPI_Finalize();
	    }
	}
	else
	{
		allocateMem(currentDataType, numElements, 0, scalars[index].data);
		

		if (currentDataType == "float")
		{
	       	MPI_File_read_all(fh, ((float *)scalars[index].data), numElements, mpiDataType, &status);

	       	float minVal, maxVal, avg;
	        minMax( ((float *)scalars[index].data), numElements, minVal, maxVal, avg);

	        minMaxAvgLog << "min: " << minVal << ", max: " << maxVal<< ", avg: " << avg << std::endl;
		}
	    else if (currentDataType == "double")
	    {
	        MPI_File_read_all(fh, ((double *)scalars[index].data), numElements, mpiDataType, &status);

	        double minVal, maxVal, avg;
	        minMax( ((double *)scalars[index].data), numElements, minVal, maxVal, avg);
	       	minMaxAvgLog << "min: " << minVal << ", max: " << maxVal<< ", avg: " << avg << std::endl;
	    }
	    else if (currentDataType == "int")
	    {
	        MPI_File_read_all(fh, ((int *)scalars[index].data), numElements, mpiDataType, &status);

	        int minVal, maxVal, avg;
	        minMax( ((int *)scalars[index].data), numElements, minVal, maxVal, avg);
	        minMaxAvgLog << "min: " << minVal << ", max: " << maxVal<< ", avg: " << avg << std::endl;
	    }
	    else
	    {	
	    	if (myRank == 0)
	    		std::cout << "Unsupported data type " << currentDataType << " . The program will now exit!" << std::endl;

	    	MPI_Finalize();
	    }
	}

	MPI_File_close( &fh );

	

	clock.stop();
	log << "\n--------------------------" << std::endl;
	log << "Param: " << paramName << std::endl;
	log << "mpiDivisions: " << mpiDivisions[0] << ", " << mpiDivisions[1] << ", " << mpiDivisions[2] << std::endl;
	log << "sizePerDim: "	<< sizePerDim[0]	<< ", " << sizePerDim[1]<< ", " << sizePerDim[2]<< std::endl;
	log << "rankOffset: "	<< rankOffset[0]	<< ", " << rankOffset[1]<< ", " << rankOffset[2]<< std::endl;
	log << minMaxAvgLog.str() << std::endl;
	log << "numElements: "	<< numElements << std::endl;
	log << "totalNumberOfElements: " << totalNumberOfElements << std::endl;

	log << "Loading data took: " << clock.getDuration() << " s" << std::endl;

	return 1;
}




inline int GenericBinaryLoader::saveCompData(std::string paramName, void * cData)
{
	Timer clock;
	clock.start();
	
	int index = 0;
	for (auto s: scalars)
	{
		if (paramName == s.name)
			break;

		index = index + 1;
	}


	if (index != scalars.size())
	{
		scalars[index].compressed = true;
		allocateMem(scalars[index].type, numElements, 0, scalars[index].data);
		memcpy(scalars[index].data, cData, numElements*getDataypeSize(scalars[index].type));
	}
	else
		return -1;

	clock.stop();
	log << "Saving data " << paramName << " took: " << clock.getDuration() << " s" << std::endl;

	return 0;
}



inline int GenericBinaryLoader::writeData(std::string _filename)
{
	Timer clock;
	clock.start();

	MPI_File fh;
	MPI_Status status;
	MPI_Datatype filetype;
	std::stringstream minMaxAvgLog;

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


	for (int i=0; i<scalars.size(); i++)
	{
		if (scalars[i].compressed == false)
		{
			scalars[i].toBeCompressed = false;
			loadData(scalars[i].name);
		}


		MPI_Offset mpi_offset = scalars[i].offset;
		MPI_File_set_view(fh, mpi_offset, MPI_BYTE, filetype, "native", MPI_INFO_NULL);

		if (scalars[i].type == "float")
		{
	       	MPI_File_write_all(fh, (float *)scalars[i].data, numElements, mpiDataType, &status);
		}
	    else if (scalars[i].type == "double")
	    {
	    	double minVal, maxVal, avg;
	        minMax( (double *)scalars[i].data, numElements, minVal, maxVal, avg);
	        minMaxAvgLog << "!!!! min: " << minVal << ", max: " << maxVal<< ", avg: " << avg << std::endl;

	        MPI_File_write_all(fh, (double *)scalars[i].data, numElements, mpiDataType, &status);  
	    }
	    else if (scalars[i].type == "int")
	    {
	        MPI_File_write_all(fh, (int *)scalars[i].data, numElements, mpiDataType, &status);
	    }

	    deAllocateMem(scalars[i].type, scalars[i].data);
	}

	MPI_File_close( &fh );

	
	clock.stop();
	log << minMaxAvgLog.str() << std::endl;
	log << "writing data took: " << clock.getDuration() << " s" << std::endl;

	return 0;
}

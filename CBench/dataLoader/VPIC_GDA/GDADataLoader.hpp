/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/

#ifndef _GDA_LOADER_H_
#define _GDA_LOADER_H_

#include <sstream>
#include <string>
#include <unordered_map>

#include "dataLoaderInterface.hpp"

#include "strConvert.hpp"
#include "json.hpp"
#include "timer.hpp"
#include "utils.hpp"


class GDADataLoader: public DataLoaderInterface
{
	int numRanks;
	int myRank;

	float origRealDims[3];

	void *tempData;
	MPI_Datatype mpiDataType;
	int mpiDivisions[3];

  public:
	GDADataLoader();
	~GDADataLoader();

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int saveCompData(std::string paramName, void * cData);
	int writeData(std::string _filename);

	int saveInputFileParameters() { return 1; };
	int close() { return deAllocateMem(dataType, data); }
	void setParam(std::string paramName, std::string type, std::string value) {};
	bool loadUncompressedFields(nlohmann::json const&) { return false; }
};


inline GDADataLoader::GDADataLoader()
{

}


inline GDADataLoader::~GDADataLoader()
{

}


inline void GDADataLoader::init(std::string _filename, MPI_Comm _comm)
{
	filename = _filename;
	comm = _comm;
	saveData = false;

	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);

	std::string metadataFile = filename + ".info";

	std::string line;
	std::ifstream myfile(metadataFile.c_str());
	if (myfile.is_open())
	{
		getline(myfile, line); // padding
		getline(myfile, line); strConvert::to_x(line, origDims[0]); // nx
		getline(myfile, line); strConvert::to_x(line, origDims[1]); // ny
		getline(myfile, line); strConvert::to_x(line, origDims[2]); // nz
		getline(myfile, line); // padding
		getline(myfile, line); strConvert::to_x(line, origRealDims[0]);// lx
		getline(myfile, line); strConvert::to_x(line, origRealDims[1]);// ly
		getline(myfile, line); strConvert::to_x(line, origRealDims[2]);// lz
		getline(myfile, dataType); 

		myfile.close();
	}
	else
	{
		if (myRank == 0)
			std::cout << "Unable to open file " << metadataFile << ". This program will now exit!" << std::endl;

		MPI_Finalize();
	}


	debugLog << "original dims: " << origDims[0] << ", " << origDims[1] << ", " << origDims[1] << std::endl;
	debugLog << "real dims : " << origRealDims[0] << ", " << origRealDims[1] << ", " << origRealDims[1] << std::endl;
	debugLog << "dataType : " << dataType << std::endl;
}



inline int GDADataLoader::loadData(std::string paramName)
{
	Timer clock("load");

	float minVal, maxVal, avg;
	totalNumberOfElements = origDims[0] * origDims[1] * origDims[2];

	//
	// TODO: MPI split!!!
	Partition current = getPartition(myRank, numRanks, origDims[0], origDims[1], origDims[2]);
	rankOffset[0] = current.min_x;
	rankOffset[1] = current.min_y;
	rankOffset[2] = current.min_z;

	sizePerDim[0] = current.max_x - current.min_x;
	sizePerDim[1] = current.max_y - current.min_y;
	sizePerDim[2] = current.max_z - current.min_z;

	numElements = sizePerDim[0] * sizePerDim[1] * sizePerDim[2];


	// Create space for data and store data there
	allocateMem(dataType, numElements, 0, data);
	elemSize = Memory::sizeOf[dataType];

	// Read the binary file
	std::string dataFile = filename + ".gda";



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


	if (dataType == "float")
       	MPI_File_read_all(fh, ((float *)data), numElements, mpiDataType, &status);
    else if (dataType == "double")
        MPI_File_read_all(fh, ((double *)data), numElements, mpiDataType, &status);
    else if (dataType == "int")
        MPI_File_read_all(fh, ((int *)data), numElements, mpiDataType, &status);
    else
    {	
    	if (myRank == 0)
    		std::cout << "Unsupported data type " << dataType << " for VPIC GDA. The program will now exit!" << std::endl;

    	MPI_Finalize();
    }

	MPI_File_close( &fh );


	// TODO: Not extensible
	minMax( ((float *)data), numElements, minVal, maxVal, avg);
	

	clock.stop("load");
	debugLog << "min: " << minVal << ", max: " << maxVal<< ", avg: " << avg << std::endl;
	debugLog << "mpiDivisions: " << mpiDivisions[0] << ", " << mpiDivisions[1] << ", " << mpiDivisions[2] << std::endl;
	debugLog << "sizePerDim: "	<< sizePerDim[0]	<< ", " << sizePerDim[1]<< ", " << sizePerDim[2]<< std::endl;
	debugLog << "rankOffset: "	<< rankOffset[0]	<< ", " << rankOffset[1]<< ", " << rankOffset[2]<< std::endl;
	debugLog << "numElements: "	<< numElements << std::endl;
	debugLog << "totalNumberOfElements: " << totalNumberOfElements << std::endl;

	debugLog << "Loading data took: " << clock.getDuration("load") << " s" << std::endl;
}



inline int GDADataLoader::saveCompData(std::string paramName, void * cData)
{
	Timer clock("save");
	
	allocateMem(dataType, numElements, 0, tempData);
	memcpy(tempData, cData, numElements*getDataypeSize(dataType));

	clock.stop("save");
	debugLog << "saving data took: " << clock.getDuration("save") << " s" << std::endl;
}



inline int GDADataLoader::writeData(std::string _filename)
{
	Timer clock("write");

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

	_filename = _filename + ".gda";
	MPI_File_open(comm, _filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	MPI_File_set_view(fh, 0, mpiDataType, filetype, "native", MPI_INFO_NULL);
	

	if (dataType == "float")
       	MPI_File_write_all(fh, (float *)tempData, numElements, mpiDataType, &status);
    else if (dataType == "double")
        MPI_File_write_all(fh, (double *)tempData, numElements, mpiDataType, &status);
    else if (dataType == "int")
        MPI_File_write_all(fh, (int *)tempData, numElements, mpiDataType, &status);
    else
    {	
    	if (myRank == 0)
    		std::cout << "Unsupported data type " << dataType << " for VPIC GDA. The program will now exit!" << std::endl;

    	MPI_Finalize();
    }

	MPI_File_close( &fh );

	deAllocateMem(dataType, tempData);

	clock.stop("write");
	debugLog << "writing data took: " << clock.getDuration("write") << " s" << std::endl;
}


#endif


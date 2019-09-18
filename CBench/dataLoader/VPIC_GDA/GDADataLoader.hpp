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
	int mpiDivisions[3];

	void *tempData;

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

		myfile.close();
	}
	else
		std::cout << "Unable to open file " << metadataFile << std::endl;

	dataType = "float";


	log << "original dims: " << origDims[0] << ", " << origDims[1] << ", " << origDims[1] << std::endl;
	log << "real dims : " << origRealDims[0] << ", " << origRealDims[1] << ", " << origRealDims[1] << std::endl;
}



inline int GDADataLoader::loadData(std::string paramName)
{
	Timer clock;
	clock.start();


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

	// Read the binary file
	std::string dataFile = filename + ".gda";
	ifstream myFile (dataFile.c_str(), ios::in | ios::binary);

	int dataSizeInBytes = 4;
	size_t indexLocal = 0;
	for (size_t z = 0; z < sizePerDim[2]; z++)
		for (size_t y = 0; y < sizePerDim[1]; y++)
		{
			size_t globalStartIndex = z * (origDims[0] * origDims[1]) * dataSizeInBytes +
			                          y * origDims[0] * dataSizeInBytes +
			                          rankOffset[0] * dataSizeInBytes;

			size_t localStartIndex = z * (sizePerDim[1] * sizePerDim[0]) +
			                         y * sizePerDim[0];

			size_t readLength = sizePerDim[2] * dataSizeInBytes;
			char * memblock = new char[readLength];

			myFile.seekg(globalStartIndex);
			myFile.read(memblock, readLength);

			memcpy ( &((( float *)data)[localStartIndex]), &memblock, readLength);

			delete [] memblock;
		}

	myFile.close();


	clock.stop();
	log << "origDims: "     << origDims[0]		<< ", " << origDims[1]<< ", " << origDims[2]<< std::endl;
	log << "sizePerDim: "	<< sizePerDim[0]	<< ", " << sizePerDim[1]<< ", " << sizePerDim[2]<< std::endl;
	log << "rankOffset: "	<< rankOffset[0]	<< ", " << rankOffset[1]<< ", " << rankOffset[2]<< std::endl;
	log << "numElements: "	<< numElements << std::endl;
	log << "totalNumberOfElements: " << totalNumberOfElements << std::endl;

	log << "Loading data took: " << clock.getDuration() << " s" << std::endl;
}



inline int GDADataLoader::saveCompData(std::string paramName, void * cData)
{
	Timer clock;
	clock.start();

	allocateMem(dataType, numElements, 0, tempData);
	memcpy(tempData, cData, numElements*getDataypeSize(dataType));

	clock.stop();
	log << "saving data took: " << clock.getDuration() << " s" << std::endl;
}



inline int GDADataLoader::writeData(std::string _filename)
{
	Timer clock;
	clock.start();

	MPI_File fh;
	MPI_Comm commWrite;
	MPI_Status status;

	int periods[3] = { 0, 0, 0 };
	int divisions[3] = {2, 2, 2};
	MPI_Cart_create(comm, 3, divisions, periods, 0, &commWrite);

	int fullsize[3], currentSize[3], offsets[3];
	for (int i = 0; i < 3; i++)
	{
		fullsize[i] = origDims[i];
		currentSize[i] = sizePerDim[i];
		offsets[i] = rankOffset[i];
	}
	MPI_Datatype filetype;
	MPI_Type_create_subarray(3, fullsize, currentSize, offsets, MPI_ORDER_C, MPI_FLOAT, &filetype);

	MPI_Type_commit(&filetype);
	MPI_File_open(commWrite, _filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	MPI_File_set_view(fh, 0, MPI_FLOAT, filetype, "native", MPI_INFO_NULL);
	MPI_File_write_all(fh, (float *)tempData, numElements, MPI_FLOAT, &status);

	deAllocateMem(dataType, tempData);

	clock.stop();
	log << "writing data took: " << clock.getDuration() << " s" << std::endl;
}


#endif


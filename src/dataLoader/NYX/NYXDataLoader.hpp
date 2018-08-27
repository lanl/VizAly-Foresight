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
	/*// Allocate mem
	if (dataType == "float")
		data = new float[numElements + offset];
	else if (dataType == "double")
		data = new double[numElements + offset];
	else if (dataType == "int8_t")
		data = new int8_t[numElements + offset];
	else if (dataType == "int16_t")
		data = new int16_t[numElements + offset];
	else if (dataType == "int32_t")
		data = new int32_t[numElements + offset];
	else if (dataType == "int64_t")
		data = new int64_t[numElements + offset];
	else if (dataType == "uint8_t")
		data = new uint8_t[numElements + offset];
	else if (dataType == "uint16_t")
		data = new uint16_t[numElements + offset];
	else if (dataType == "uint32_t")
		data = new uint32_t[numElements + offset];
	else if (dataType == "uint64_t")
		data = new uint64_t[numElements + offset];
	else
		return 0;

	// Get size
	if (dataType == "float")
		elemSize = sizeof(float);
	else if (dataType == "double")
		elemSize = sizeof(double);
	else if (dataType == "int8_t")
		elemSize = sizeof(int8_t);
	else if (dataType == "int16_t")
		elemSize = sizeof(int16_t);
	else if (dataType == "int32_t")
		elemSize = sizeof(int32_t);
	else if (dataType == "int64_t")
		elemSize = sizeof(int64_t);
	else if (dataType == "uint8_t")
		elemSize = sizeof(uint8_t);
	else if (dataType == "uint16_t")
		elemSize = sizeof(uint16_t);
	else if (dataType == "uint32_t")
		elemSize = sizeof(uint32_t);
	else if (dataType == "uint64_t")
		elemSize = sizeof(uint64_t);
	else
		return 0;*/

	return 1;
}


inline int NYXDataLoader::deAllocateMem()
{
	if (data == NULL) // already deallocated!
		return 1;
	/*
	if (dataType == "float")
		delete[](float*) data;
	else if (dataType == "double")
		delete[](double*) data;
	else if (dataType == "int8_t")
		delete[](int8_t*) data;
	else if (dataType == "int16_t")
		delete[](int16_t*) data;
	else if (dataType == "int32_t")
		delete[](int32_t*) data;
	else if (dataType == "int64_t")
		delete[](int64_t*) data;
	else if (dataType == "uint8_t")
		delete[](uint8_t*) data;
	else if (dataType == "uint16_t")
		delete[](uint16_t*) data;
	else if (dataType == "uint32_t")
		delete[](uint32_t*) data;
	else if (dataType == "uint64_t")
		delete[](uint64_t*) data;
	else
		return 0;*/

	data = NULL;

	return 1;
}


inline int NYXDataLoader::loadData(std::string paramName)
{
	Timer clock;
	log.str("");

	clock.start();
	/*gio::GenericIO *gioReader;
	param = paramName;

	// Init GenericIO reader + open file
	gioReader = new gio::GenericIO(filename, gio::GenericIO::FileIOPOSIX);

	// Open file
	gioReader->openAndReadHeader(gio::GenericIO::MismatchRedistribute);
	int numDataRanks = gioReader->readNRanks();

	if (numRanks > numDataRanks)
	{
		std::cout << "Num data ranks: " << numDataRanks << "Use <= MPI ranks than data ranks" << std::endl;
		return -1;
	}

	// Count number of elements
	totalNumberOfElements = 0;
	for (int i = 0; i < numDataRanks; ++i)
		totalNumberOfElements += gioReader->readNumElems(i);

	// Read in the scalars information
	std::vector<gio::GenericIO::VariableInfo> VI;
	gioReader->getVariableInfo(VI);
	int numVars = static_cast<int>(VI.size());

	// Store scalar info
	std::vector<GioData> readInData;
	readInData.resize(numVars);

	int paramToLoadPos = -1;
	for (int i = 0; i < numVars; i++)
	{
		readInData[i].id = i;
		readInData[i].name = VI[i].Name;
		readInData[i].size = static_cast<int>(VI[i].Size);
		readInData[i].isFloat = VI[i].IsFloat;
		readInData[i].isSigned = VI[i].IsSigned;
		readInData[i].ghost = VI[i].MaybePhysGhost;
		readInData[i].xVar = VI[i].IsPhysCoordX;
		readInData[i].yVar = VI[i].IsPhysCoordY;
		readInData[i].zVar = VI[i].IsPhysCoordZ;
		readInData[i].determineDataType();

		// Check if this is one of the scalars to read
		if (readInData[i].name == paramName)
		{
			readInData[i].loadData = true;
			paramToLoadPos = i;
		}
	}

	if (paramToLoadPos == -1)
	{
		std::cout << "Cannot find that parameter, exiting now!";
		return -2;
	}


	//
	// Split ranks among data
	int numDataRanksPerMPIRank = numDataRanks / numRanks;
	int loadRange[2];
	loadRange[0] = myRank * numDataRanksPerMPIRank;
	loadRange[1] = (myRank + 1) * numDataRanksPerMPIRank;
	if (myRank == numRanks - 1)
		loadRange[1] = numDataRanks;


	//
	// Determine memory size and allocate memory
	readInData[paramToLoadPos].determineDataType();
	dataType = readInData[paramToLoadPos].dataType;

	numElements = 0;
	for (int i = loadRange[0]; i < loadRange[1]; i++)
		numElements += gioReader->readNumElems(i);

	allocateMem(dataType, numElements, 0);


	//
	// Actually load the data
	size_t offset = 0;
	size_t offsetSize = readInData[paramToLoadPos].size;


	for (int i = loadRange[0]; i < loadRange[1]; i++) // for each rank
	{
		size_t Np = gioReader->readNumElems(i);

		if (readInData[paramToLoadPos].dataType == "float")
			gioReader->addVariable( (readInData[paramToLoadPos].name).c_str(), (float*)data + offset, true);
		else if (readInData[paramToLoadPos].dataType == "double")
			gioReader->addVariable( (readInData[paramToLoadPos].name).c_str(), (double*)data + offset, true);
		else if (readInData[paramToLoadPos].dataType == "int8_t")
			gioReader->addVariable((readInData[paramToLoadPos].name).c_str(), (int8_t*)data + offset, true);
		else if (readInData[paramToLoadPos].dataType == "int16_t")
			gioReader->addVariable((readInData[paramToLoadPos].name).c_str(), (int16_t*)data + offset, true);
		else if (readInData[paramToLoadPos].dataType == "int32_t")
			gioReader->addVariable((readInData[paramToLoadPos].name).c_str(), (int32_t*)data + offset, true);
		else if (readInData[paramToLoadPos].dataType == "int64_t")
			gioReader->addVariable((readInData[paramToLoadPos].name).c_str(), (int64_t*)data + offset, true);
		else if (readInData[paramToLoadPos].dataType == "uint8_t")
			gioReader->addVariable((readInData[paramToLoadPos].name).c_str(), (uint8_t*)data + offset, true);
		else if (readInData[paramToLoadPos].dataType == "uint16_t")
			gioReader->addVariable((readInData[paramToLoadPos].name).c_str(), (uint16_t*)data + offset, true);
		else if (readInData[paramToLoadPos].dataType == "uint32_t")
			gioReader->addVariable((readInData[paramToLoadPos].name).c_str(), (uint32_t*)data + offset, true);
		else if (readInData[paramToLoadPos].dataType == "uint64_t")
			gioReader->addVariable((readInData[paramToLoadPos].name).c_str(), (uint64_t*)data + offset, true);

		else
			std::cout << " = data type undefined!!!" << std::endl;

		gioReader->readDataSection(0, Np, i, false); // reading the whole file

		offset = offset + Np;
	}
	clock.stop();*/
	log << "Loading data took " << clock.getDuration() << " s" << std::endl;

	return 1; // All good
}

#endif
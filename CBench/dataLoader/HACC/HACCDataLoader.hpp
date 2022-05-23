/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/
#pragma once

#include <algorithm>
#include <sstream>
#include <climits>
#include <memory>
#include "thirdparty/genericio/GenericIO.h"
#include "gioData.hpp"
#include "dataLoaderInterface.hpp"
#include "timer.hpp"
#include "utils.hpp"


class HACCDataLoader: public DataLoaderInterface
{
	int numRanks;
	int myRank;

	// For output
	double physOrigin[3];
	double physScale[3];
	int mpiCartPartitions[3];

  public:
	HACCDataLoader();
	~HACCDataLoader();

	int loadData(std::string paramName, void *& _data);

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int saveCompData(std::string paramName, void * cData);
	int writeData(std::string _filename);
	int saveInputFileParameters();
	int close() { return deAllocateMem(dataType, data); }
	void setParam(std::string paramName, std::string type, std::string value){};
  bool loadUncompressedFields(nlohmann::json const&) { return false; }
};


inline HACCDataLoader::HACCDataLoader()
{
	myRank = 0;
	numRanks = 0;
	loader = "HACC";
	saveData = false;

	inOutData.clear();
}


inline HACCDataLoader::~HACCDataLoader()
{
	deAllocateMem(dataType, data);
}


inline void HACCDataLoader::init(std::string _filename, MPI_Comm _comm)
{
	filename = _filename;
	comm = _comm;
	saveData = false;

	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);
}




inline int HACCDataLoader::saveInputFileParameters()
{
	std::unique_ptr<gio::GenericIO> gioReader( new gio::GenericIO(comm, filename));

	// Open filetotalNumberOfElements =
	gioReader->openAndReadHeader(gio::GenericIO::MismatchRedistribute);
	
	int numDataRanks = gioReader->readNRanks();

	if (numRanks > numDataRanks)
	{
		std::cout << "Num data ranks: " << numDataRanks << "Use <= MPI ranks than data ranks" << std::endl;
		return -1;
	}


	// Get dimensions of the input file
	gioReader->readPhysOrigin(physOrigin);
    gioReader->readPhysScale(physScale);


	// Read in the scalars information
	std::vector<gio::GenericIO::VariableInfo> VI;
	gioReader->getVariableInfo(VI);
	int numVars = static_cast<int>(VI.size());

	
	for (int i = 0; i < numVars; i++)
	{
		GioData readInData(i, VI[i].Name, static_cast<int>(VI[i].Size), VI[i].IsFloat, VI[i].IsSigned, VI[i].IsPhysCoordX, VI[i].IsPhysCoordY, VI[i].IsPhysCoordZ);
		readInData.determineDataType();

		inOutData.push_back(readInData);
	}

	return 1;
}


inline int HACCDataLoader::loadData(std::string paramName)
{
	Timer clock;
	

	clock.start("load");
	gio::GenericIO *gioReader;
	param = paramName;

	// Init GenericIO reader + open file
	gioReader = new gio::GenericIO(comm, filename);
	//gioReader = new gio::GenericIO(filename);




	// Open file
	gioReader->openAndReadHeader(gio::GenericIO::MismatchRedistribute);
	int numDataRanks = gioReader->readNRanks();

	if (numRanks > numDataRanks)
	{
		std::cout << "Num data ranks: " << numDataRanks << "Use <= MPI ranks than data ranks" << std::endl;
		return -1;
	}

	debugLog << "numRanks: " << numRanks << ", numDataRanks: " << numDataRanks << std::endl;

	// Count number of elements
	totalNumberOfElements = 0;
	for (int i = 0; i < numDataRanks; ++i)
		totalNumberOfElements += gioReader->readNumElems(i);


	// Read in the scalars information
	std::vector<gio::GenericIO::VariableInfo> VI;
	gioReader->getVariableInfo(VI);
	int numVars = static_cast<int>(VI.size());


	bool paramToLoad = false;
	GioData readInData;
	for (int i = 0; i < numVars; i++)
		if (VI[i].Name == paramName)
		{
			paramToLoad = true;
			readInData.init(i, VI[i].Name, static_cast<int>(VI[i].Size), VI[i].IsFloat, VI[i].IsSigned, VI[i].IsPhysCoordX, VI[i].IsPhysCoordY, VI[i].IsPhysCoordZ);
			readInData.determineDataType();
		
			dataType = readInData.dataType;
			elemSize = readInData.size;

			break;
		}


	if ( !paramToLoad )
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

	int splitDims[3];
	gioReader->readDims(splitDims);
	debugLog << "splitDims: " << splitDims[0] << "," << splitDims[1] << "," << splitDims[2] << std::endl;
	debugLog << myRank << " ~ loadRange[0]: " << loadRange[0] << ", loadRange[1]: " << loadRange[1] << std::endl;







	//
	// Determine memory size and allocate memory
	size_t maxNumElementsPerRank = 0;
	numElements = 0;
	for (int i = loadRange[0]; i < loadRange[1]; i++)
	{
		numElements += gioReader->readNumElems(i);
		maxNumElementsPerRank = std::max(maxNumElementsPerRank,numElements);
	}


	allocateMem(dataType, numElements, 0, data);


	readInData.setNumElements(maxNumElementsPerRank);
	readInData.alloc(1);

	sizePerDim[0] = numElements;	// For compression



	
	//
	// Actually load the data
	int minmaxX[2] = {INT_MAX, INT_MIN};
	int minmaxY[2] = {INT_MAX, INT_MIN};
	int minmaxZ[2] = {INT_MAX, INT_MIN};

	size_t offset = 0;
	for (int i = loadRange[0]; i < loadRange[1]; i++) // for each rank
	{
		size_t Np = gioReader->readNumElems(i);

		int coords[3];
		gioReader->readCoords(coords, i);
		debugLog << "Coord indices: " << coords[0] << ", " << coords[1] << ", " << coords[2] << " | ";

		debugLog << "coordinates: (" << (float)coords[0] / splitDims[0] * physScale[0] + physOrigin[0] << ", "
                          << (float)coords[1] / splitDims[1] * physScale[1] + physOrigin[1] << ", "
                          << (float)coords[2] / splitDims[2] * physScale[2] + physOrigin[2] << ") -> ("
                          << (float)(coords[0] + 1) / splitDims[0] * physScale[0] + physOrigin[0] << ", "
                          << (float)(coords[1] + 1) / splitDims[1] * physScale[1] + physOrigin[1] << ", "
                          << (float)(coords[2] + 1) / splitDims[2] * physScale[2] + physOrigin[2] << ")" << std::endl;

        if (saveData)
        {
        	minmaxX[0] = std::min( (int) ((float)coords[0] / splitDims[0] * physScale[0] + physOrigin[0]), minmaxX[0] );
        	minmaxY[0] = std::min( (int) ((float)coords[1] / splitDims[1] * physScale[1] + physOrigin[1]), minmaxY[0] );
        	minmaxZ[0] = std::min( (int) ((float)coords[2] / splitDims[2] * physScale[2] + physOrigin[2]), minmaxZ[0] );

        	minmaxX[1] = std::max( (int) ((float)(coords[0] + 1) / splitDims[0] * physScale[0] + physOrigin[0]), minmaxX[1] );
        	minmaxY[1] = std::max( (int) ((float)(coords[1] + 1) / splitDims[1] * physScale[1] + physOrigin[1]), minmaxY[1] );
        	minmaxZ[1] = std::max( (int) ((float)(coords[2] + 1) / splitDims[2] * physScale[2] + physOrigin[2]), minmaxZ[1] );
        }
		
	
		if (readInData.dataType == "float")
			gioReader->addVariable( (readInData.name).c_str(), (float*)readInData.data, true);
		else if (readInData.dataType == "double")
			gioReader->addVariable( (readInData.name).c_str(), (double*)readInData.data, true);
		else if (readInData.dataType == "int8_t")
			gioReader->addVariable((readInData.name).c_str(), (int8_t*)readInData.data, true);
		else if (readInData.dataType == "int16_t")
			gioReader->addVariable((readInData.name).c_str(), (int16_t*)readInData.data, true);
		else if (readInData.dataType == "int32_t")
			gioReader->addVariable((readInData.name).c_str(), (int32_t*)readInData.data, true);
		else if (readInData.dataType == "int64_t")
			gioReader->addVariable((readInData.name).c_str(), (int64_t*)readInData.data, true);
		else if (readInData.dataType == "uint8_t")
			gioReader->addVariable((readInData.name).c_str(), (uint8_t*)readInData.data, true);
		else if (readInData.dataType == "uint16_t")
			gioReader->addVariable((readInData.name).c_str(), (uint16_t*)readInData.data, true);
		else if (readInData.dataType == "uint32_t")
			gioReader->addVariable((readInData.name).c_str(), (uint32_t*)readInData.data, true);
		else if (readInData.dataType == "uint64_t")
			gioReader->addVariable((readInData.name).c_str(), (uint64_t*)readInData.data, true);
	
		//gioReader->readDataSection(0, Np, i, false); // reading the whole file
		gioReader->readData(i, false); // reading the whole file


		if (readInData.dataType == "float")
			memcpy( &((float*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "double")
			memcpy( &((double*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "int8_t")
			memcpy( &((int8_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "int16_t")
			memcpy( &((int16_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "int32_t")
			memcpy( &((int32_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "int64_t")
			memcpy( &((int64_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "uint8_t")
			memcpy( &((uint8_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "uint16_t")
			memcpy( &((uint16_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "uint32_t")
			memcpy( &((uint32_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "uint64_t")
			memcpy( &((uint64_t*)data)[offset],  readInData.data, Np*readInData.size);
		
		
		offset = offset + Np;
	}
	clock.stop("load");


	if (saveData)
	{
		int rangeX = minmaxX[1]-minmaxX[0];
		int rangeY = minmaxY[1]-minmaxY[0];
		int rangeZ = minmaxZ[1]-minmaxZ[0];

		mpiCartPartitions[0] = physScale[0]/rangeX;
		mpiCartPartitions[1] = physScale[1]/rangeY;
		mpiCartPartitions[2] = physScale[2]/rangeZ;

		debugLog  << "mpiCartPartitions: " << mpiCartPartitions[0] << ", " << mpiCartPartitions[1] << ", " << mpiCartPartitions[2] << std::endl;
	}
	
	deAllocateMem(dataType, readInData.data);


	return 1; // All good
}

/*
inline int HACCDataLoader::loadData(std::string paramName)
{
	Timer clock;
	clock.start("load");

	param = paramName;

	// Init GenericIO reader + open file
	std::unique_ptr<gio::GenericIO> gioReader( new gio::GenericIO(comm, filename));

	gio::GenericIO::setNaturalDefaultPartition();

	// Open file
	gioReader->openAndReadHeader(gio::GenericIO::MismatchRedistribute);
	int numDataRanks = gioReader->readNRanks();
	totalNumberOfElements = gioReader->readTotalNumElems();

	if (numRanks > numDataRanks)
	{
		std::cout << "Num data ranks: " << numDataRanks << "Use <= MPI ranks than data ranks" << std::endl;
		return -1;
	}

	std::cout << myRank << " ~ num MPI Ranks: " << numRanks << ", num Data Ranks: " << numDataRanks << std::endl;
	debugLog << "num MPI Ranks: " << numRanks << ", num Data Ranks: " << numDataRanks << std::endl;

	//
	// Split ranks among data
	int numDataRanksPerMPIRank = numDataRanks / numRanks;
	int loadRange[2];
	loadRange[0] = myRank * numDataRanksPerMPIRank;
	loadRange[1] = (myRank + 1) * numDataRanksPerMPIRank;
	if (myRank == numRanks - 1)
		loadRange[1] = numDataRanks;

	debugLog << myRank << " ~ loadRange[0]: " << loadRange[0] << ", loadRange[1]: " << loadRange[1] << std::endl;


	// MPI_Barrier(comm);
	// std::cout << myRank << " ~ loadRange[0]: " << loadRange[0] << ", loadRange[1]: " << loadRange[1] << std::endl;
	// MPI_Barrier(comm);


	//
	// Determine memory size and allocate memory
	numElements = 0;
	for (int i = loadRange[0]; i < loadRange[1]; i++)
		numElements += gioReader->readNumElems(i);

	// MPI_Barrier(comm);
	// std::cout << myRank << " ~ num Elements in rank: " << numElements << std::endl;
	// MPI_Barrier(comm);



	// Read in the scalars information
	std::vector<gio::GenericIO::VariableInfo> VI;
	gioReader->getVariableInfo(VI);
	int numVars = static_cast<int>(VI.size());


	bool paramToLoad = false;
	GioData readInData;
	for (int i = 0; i < numVars; i++)
		if (VI[i].Name == paramName)
		{
			paramToLoad = true;
			readInData.init(i, VI[i].Name, static_cast<int>(VI[i].Size), VI[i].IsFloat, VI[i].IsSigned, VI[i].IsPhysCoordX, VI[i].IsPhysCoordY, VI[i].IsPhysCoordZ);
			readInData.determineDataType();
		
			dataType = readInData.dataType;
			elemSize = readInData.size;

			break;
		}


	if ( !paramToLoad )
	{
		std::cout << "Cannot find that parameter, exiting now!";
		return -2;
	}


	allocateMem(dataType, numElements, 0, data);


	readInData.setNumElements(numElements);
	readInData.alloc(1);

	sizePerDim[0] = numElements;	// For compression

	std::cout << myRank << " ~ loading " << numDataRanksPerMPIRank << " data ranks ... "<< std::endl;

	//
	// Actually load the data
	int minmaxX[2] = {INT_MAX, INT_MIN};
	int minmaxY[2] = {INT_MAX, INT_MIN};
	int minmaxZ[2] = {INT_MAX, INT_MIN};

	int splitDims[3];
	gioReader->readDims(splitDims);

	
	debugLog << "splitDims: " << splitDims[0] << "," << splitDims[1] << "," << splitDims[2] << std::endl;
	

	size_t offset = 0;
	for (int i = loadRange[0]; i < loadRange[1]; i++) // for each rank
	{
		size_t Np = gioReader->readNumElems(i);

		int coords[3];
		gioReader->readCoords(coords, i);
		debugLog << "Coord indices: " << coords[0] << ", " << coords[1] << ", " << coords[2] << " | ";

		debugLog << "coordinates: (" << (float)coords[0] / splitDims[0] * physScale[0] + physOrigin[0] << ", "
                          << (float)coords[1] / splitDims[1] * physScale[1] + physOrigin[1] << ", "
                          << (float)coords[2] / splitDims[2] * physScale[2] + physOrigin[2] << ") -> ("
                          << (float)(coords[0] + 1) / splitDims[0] * physScale[0] + physOrigin[0] << ", "
                          << (float)(coords[1] + 1) / splitDims[1] * physScale[1] + physOrigin[1] << ", "
                          << (float)(coords[2] + 1) / splitDims[2] * physScale[2] + physOrigin[2] << ")" << std::endl;

        if (saveData)
        {
        	minmaxX[0] = std::min( (int) ((float)coords[0] / splitDims[0] * physScale[0] + physOrigin[0]), minmaxX[0] );
        	minmaxY[0] = std::min( (int) ((float)coords[1] / splitDims[1] * physScale[1] + physOrigin[1]), minmaxY[0] );
        	minmaxZ[0] = std::min( (int) ((float)coords[2] / splitDims[2] * physScale[2] + physOrigin[2]), minmaxZ[0] );

        	minmaxX[1] = std::max( (int) ((float)(coords[0] + 1) / splitDims[0] * physScale[0] + physOrigin[0]), minmaxX[1] );
        	minmaxY[1] = std::max( (int) ((float)(coords[1] + 1) / splitDims[1] * physScale[1] + physOrigin[1]), minmaxY[1] );
        	minmaxZ[1] = std::max( (int) ((float)(coords[2] + 1) / splitDims[2] * physScale[2] + physOrigin[2]), minmaxZ[1] );
        }
		
	
		if (readInData.dataType == "float")
			gioReader->addVariable( (readInData.name).c_str(), (float*)readInData.data, true);
		else if (readInData.dataType == "double")
			gioReader->addVariable( (readInData.name).c_str(), (double*)readInData.data, true);
		else if (readInData.dataType == "int8_t")
			gioReader->addVariable((readInData.name).c_str(), (int8_t*)readInData.data, true);
		else if (readInData.dataType == "int16_t")
			gioReader->addVariable((readInData.name).c_str(), (int16_t*)readInData.data, true);
		else if (readInData.dataType == "int32_t")
			gioReader->addVariable((readInData.name).c_str(), (int32_t*)readInData.data, true);
		else if (readInData.dataType == "int64_t")
			gioReader->addVariable((readInData.name).c_str(), (int64_t*)readInData.data, true);
		else if (readInData.dataType == "uint8_t")
			gioReader->addVariable((readInData.name).c_str(), (uint8_t*)readInData.data, true);
		else if (readInData.dataType == "uint16_t")
			gioReader->addVariable((readInData.name).c_str(), (uint16_t*)readInData.data, true);
		else if (readInData.dataType == "uint32_t")
			gioReader->addVariable((readInData.name).c_str(), (uint32_t*)readInData.data, true);
		else if (readInData.dataType == "uint64_t")
			gioReader->addVariable((readInData.name).c_str(), (uint64_t*)readInData.data, true);
	
		//gioReader->readDataSection(0, Np, i, false); // reading the whole file
		gioReader->readData(i, false); // reading the whole file


		if (readInData.dataType == "float")
			memcpy( &((float*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "double")
			memcpy( &((double*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "int8_t")
			memcpy( &((int8_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "int16_t")
			memcpy( &((int16_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "int32_t")
			memcpy( &((int32_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "int64_t")
			memcpy( &((int64_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "uint8_t")
			memcpy( &((uint8_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "uint16_t")
			memcpy( &((uint16_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "uint32_t")
			memcpy( &((uint32_t*)data)[offset],  readInData.data, Np*readInData.size);
		else if (readInData.dataType == "uint64_t")
			memcpy( &((uint64_t*)data)[offset],  readInData.data, Np*readInData.size);
		
		
		offset = offset + Np;
	}
	clock.stop("load");
	debugLog << "HACC loading data took: " << clock.getDuration("load") << " s" << std::endl;


	if (saveData)
	{
		int rangeX = minmaxX[1]-minmaxX[0];
		int rangeY = minmaxY[1]-minmaxY[0];
		int rangeZ = minmaxZ[1]-minmaxZ[0];

		mpiCartPartitions[0] = physScale[0]/rangeX;
		mpiCartPartitions[1] = physScale[1]/rangeY;
		mpiCartPartitions[2] = physScale[2]/rangeZ;

		debugLog  << "mpiCartPartitions: " << mpiCartPartitions[0] << ", " << mpiCartPartitions[1] << ", " << mpiCartPartitions[2] << std::endl;
	}
	
	deAllocateMem(dataType, readInData.data);


	return 1; // All good
}
*/


inline int HACCDataLoader::saveCompData(std::string paramName, void * cData)
{

	for (int i=0; i<inOutData.size(); i++)
	{
		if (inOutData[i].name == paramName)
		{
			inOutData[i].setNumElements(numElements);
			inOutData[i].alloc();
			memcpy(inOutData[i].data, cData, inOutData[i].size*numElements);

			inOutData[i].doWrite = true;

			//log.str("");
			debugLog << "\nHACCDataLoader::saveCompData" << std::endl;
			debugLog << paramName << " found. It has " << inOutData[i].numElements << " elements of size " << inOutData[i].size << std::endl;
		}
	}

	return 1;
}


inline int HACCDataLoader::writeData(std::string _filename)
{
  #ifndef GENERICIO_NO_MPI
	Timer clock;
	clock.start("write");
	
	// Create setup
	int periods[3] = { 0, 0, 0 };	
	MPI_Cart_create(comm, 3, mpiCartPartitions, periods, 0, &comm);


	// Init GenericIO writer + open file
	std::unique_ptr<gio::GenericIO> gioWriter( new gio::GenericIO(comm, _filename));
	gioWriter->setNumElems(numElements);


	// Init physical parameters
	for (int d = 0; d < 3; ++d)
	{
		gioWriter->setPhysOrigin(physOrigin[d], d);
		gioWriter->setPhysScale(physScale[d], d);
	}

	MPI_Barrier(comm);
	// Populate parameters
	for (int i=0; i<inOutData.size(); i++)
	{
		unsigned flag = gio::GenericIO::VarHasExtraSpace;
		if (inOutData[i].xVar)
			flag |= gio::GenericIO::VarIsPhysCoordX;
		else if (inOutData[i].yVar)
			flag |= gio::GenericIO::VarIsPhysCoordY;
		else if (inOutData[i].zVar)
			flag |= gio::GenericIO::VarIsPhysCoordZ;


		if (inOutData[i].dataType == "float")
          	gioWriter->addVariable( (inOutData[i].name).c_str(), (float*)inOutData[i].data,    flag );
        else if (inOutData[i].dataType == "double")
          	gioWriter->addVariable( (inOutData[i].name).c_str(), (double*)inOutData[i].data,   flag );
        else if (inOutData[i].dataType == "int8_t")
          	gioWriter->addVariable(  (inOutData[i].name).c_str(), (int8_t*)inOutData[i].data,   flag);
        else if (inOutData[i].dataType == "int16_t")
          	gioWriter->addVariable( (inOutData[i].name).c_str(), (int16_t*)inOutData[i].data,   flag);
        else if (inOutData[i].dataType == "int32_t")
          	gioWriter->addVariable( (inOutData[i].name).c_str(), (int32_t*)inOutData[i].data,   flag);
        else if (inOutData[i].dataType == "int64_t")
          	gioWriter->addVariable( (inOutData[i].name).c_str(), (int64_t*)inOutData[i].data,   flag);
        else if (inOutData[i].dataType == "uint8_t")
          	gioWriter->addVariable( (inOutData[i].name).c_str(), (uint8_t*)inOutData[i].data,   flag);
        else if (inOutData[i].dataType == "uint16_t")
          	gioWriter->addVariable(  (inOutData[i].name).c_str(), (uint16_t*)inOutData[i].data, flag);
        else if (inOutData[i].dataType == "uint32_t")
          	gioWriter->addVariable( (inOutData[i].name).c_str(), (uint32_t*)inOutData[i].data,  flag);
        else if (inOutData[i].dataType == "uint64_t")
          	gioWriter->addVariable( (inOutData[i].name).c_str(), (uint64_t*)inOutData[i].data,  flag);
        else
          	std::cout << " = data type undefined!!!" << std::endl;		
	}

	gioWriter->write();
	debugLog << "HACCDataLoader::writeData " << _filename << std::endl;

	MPI_Barrier(comm);

	clock.stop("write");
	debugLog << "Writing data took " << clock.getDuration("write") << " s" << std::endl;
  #endif
	return 1;
}

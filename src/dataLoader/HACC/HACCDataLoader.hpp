/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/

#ifndef _HACC_LOADER_H_
#define _HACC_LOADER_H_

#include <sstream>
#include "thirdparty/genericio/GenericIO.h"
#include "gioData.hpp"
#include "dataLoaderInterface.hpp"
#include "timer.hpp"

class HACCDataLoader: public DataLoaderInterface
{
	int numRanks;
	int myRank;

	std::vector<GioData> inOutData;

  public:
	HACCDataLoader();
	~HACCDataLoader();
	int allocateMem(std::string dataType, size_t numElements, int offset);
	int deAllocateMem();
	

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int saveCompData(std::string paramName, void * cData);
	int writeData(std::string _filename);
	int saveInputFileParameters();
	int close() { return deAllocateMem(); }
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
	deAllocateMem();

	// If write() wasnt called to flush compFullData, free it here
	for (std::pair<std::string, void *> element : compFullData)
	{
		if(element.second != NULL)
			std::free(element.second);
	}
}


inline void HACCDataLoader::init(std::string _filename, MPI_Comm _comm)
{
	filename = _filename;
	comm = _comm;
	saveData = false;

	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);
}



inline int HACCDataLoader::allocateMem(std::string dataType, size_t numElements, int offset)
{
	// Allocate mem
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
		return 0;

	return 1;
}


inline int HACCDataLoader::deAllocateMem()
{
	if (data == NULL) // already deallocated!
		return 1;

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
		return 0;

	data = NULL;

	return 1;
}



inline int HACCDataLoader::saveInputFileParameters()
{
	gio::GenericIO *gioReader;
  #ifndef GENERICIO_NO_MPI
	gioReader = new gio::GenericIO(comm, filename);
  #else
	gioReader = new gio::GenericIO(filename, gio::GenericIO::FileIOPOSIX);
  #endif

	// Open file
	gioReader->openAndReadHeader(gio::GenericIO::MismatchRedistribute);
	int numDataRanks = gioReader->readNRanks();

	if (numRanks > numDataRanks)
	{
		std::cout << "Num data ranks: " << numDataRanks << "Use <= MPI ranks than data ranks" << std::endl;
		return -1;
	}

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
	log.str("");

	clock.start();
	gio::GenericIO *gioReader;
	param = paramName;

	// Init GenericIO reader + open file
  #ifndef GENERICIO_NO_MPI
	gioReader = new gio::GenericIO(comm, filename);
  #else
	gioReader = new gio::GenericIO(filename, gio::GenericIO::FileIOPOSIX);
  #endif


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


	bool paramToLoad = false;
	GioData readInData;
	for (int i = 0; i < numVars; i++)
	{
		if (VI[i].Name == paramName)
			paramToLoad = true;
		else
			continue;

		readInData.init(i, VI[i].Name, static_cast<int>(VI[i].Size), VI[i].IsFloat, VI[i].IsSigned, VI[i].IsPhysCoordX, VI[i].IsPhysCoordY, VI[i].IsPhysCoordZ);
		readInData.determineDataType();
		dataType = readInData.dataType;
	}


	if (!paramToLoad)
	{
		std::cout << "Cannot find that parameter, exiting now!";
		return -2;
	}


	//
	// Split ranks among data
	int splitDims[3];
	gioReader->readDims(splitDims);

	int numDataRanksPerMPIRank = numDataRanks / numRanks;
	int loadRange[2];
	loadRange[0] = myRank * numDataRanksPerMPIRank;
	loadRange[1] = (myRank + 1) * numDataRanksPerMPIRank;
	if (myRank == numRanks - 1)
		loadRange[1] = numDataRanks;



	//
	// Determine memory size and allocate memory
	size_t maxNumElementsPerRank = 0;
	numElements = 0;
	for (int i = loadRange[0]; i < loadRange[1]; i++)
	{
		numElements += gioReader->readNumElems(i);
		maxNumElementsPerRank = std::max(maxNumElementsPerRank,numElements);
	}


	allocateMem(dataType, numElements, 0);

	readInData.setNumElements(maxNumElementsPerRank);
	readInData.allocateMem(1);


	// WHY ???????
	dims[0] = numElements;
	std::cout << "totalNumberOfElements" << totalNumberOfElements << std::endl;
	std::cout << "numElements" << numElements << std::endl;

	
	//
	// Actually load the data
	size_t offset = 0;
	for (int i = loadRange[0]; i < loadRange[1]; i++) // for each rank
	{
		size_t Np = gioReader->readNumElems(i);

	
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
	
		gioReader->readDataSection(0, Np, i, false); // reading the whole file

		

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
	clock.stop();

	readInData.deAllocateMem();

	std::cout << "HACCDataLoader::loadData " << paramName << " done!"  << std::endl;

	return 1; // All good
}



inline int HACCDataLoader::saveCompData(std::string paramName, void * cData)
{
	std::cout << "HACCDataLoader::saveCompData " << paramName << std::endl;
	for (int i=0; i<inOutData.size(); i++)
	{
		if (inOutData[i].name == paramName)
		{
			inOutData[i].setNumElements(totalNumberOfElements);
			inOutData[i].allocateMem();
			memcpy(inOutData[i].data, cData, 4*totalNumberOfElements);

			inOutData[i].doWrite = true;
		}
	}

	return 1;
}


inline int HACCDataLoader::writeData(std::string _filename)
{
	std::cout << "HACCDataLoader::writeData " << _filename << std::endl;
	Timer clock;
	log.str("");

	gio::GenericIO *gioWriter;

	int dims[3]= {2, 2, 2};		// 8
	int periods[3] = { 0, 0, 0 };
	MPI_Cart_create(comm, 3, dims, periods, 0, &comm);

	// Init GenericIO writer + open file
  #ifndef GENERICIO_NO_MPI
	gioWriter = new gio::GenericIO(comm, _filename);// , gio::GenericIO::FileIOMPI);
  #else
	gioWriter = new gio::GenericIO(_filename, gio::GenericIO::FileIOPOSIX);
  #endif

	gioWriter->setNumElems(numElements);



	// Init physical parameters
	int physOrigin[3] = { 0,    0,    0 };
	int physScale[3]  = { 256, 256, 256 };
	for (int d = 0; d < 3; ++d)
	{
		gioWriter->setPhysOrigin(physOrigin[d], d);
		gioWriter->setPhysScale(physScale[d], d);
	}

	// Populate parameters
	//for (std::pair<std::string, void *> element : compFullData)
	for (int i=0; i<inOutData.size(); i++)
	{
		if (inOutData[i].doWrite)
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
	}

  #ifndef GENERICIO_NO_MPI
	gioWriter->write();
	std::cout << "HACCDataLoader::writeData " << _filename << "  gioWriter->write() " << std::endl;
	MPI_Barrier(comm);
  #endif


	clock.stop();
	log << "Writing data took " << clock.getDuration() << " s" << std::endl;
	return 1;
}

#endif
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
#include <string>
#include "dataLoaderInterface.hpp"
#include "timer.hpp"
#include "H5Cpp.h"

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
	int saveCompData(std::string paramName, void * cData);
	int writeData(std::string _filename);
	int saveInputFileParameters();
	int close() { return deAllocateMem(); }
};


inline NYXDataLoader::NYXDataLoader()
{
	myRank = 0;
	numRanks = 0;
	loader = "NYX";
	saveData = false;
}

inline NYXDataLoader::~NYXDataLoader()
{
	deAllocateMem();
}


inline void NYXDataLoader::init(std::string _filename, MPI_Comm _comm)
{
	filename = _filename;
	comm = _comm;
	saveData = false;

	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);
}

inline int NYXDataLoader::saveInputFileParameters()
{
    // 

    return 1;
}

inline int NYXDataLoader::allocateMem(std::string dataType, size_t numElements, int offset)
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


inline int NYXDataLoader::deAllocateMem()
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


inline int NYXDataLoader::loadData(std::string paramName)
{
	Timer clock;
	log.str("");

	clock.start();
	
	// Note: hdf5-cxx not compatible with parallel
	try {
		H5::H5File file(filename, H5F_ACC_RDONLY);
		H5::Group group(file.openGroup("native_fields"));
		/*for (int grps = 0; grps < group.getNumObjs(); grps++)
		{
			std::cout << "Field: " << group.getObjnameByIdx(grps) << "\n";
			std::string name = group.getObjnameByIdx(grps);
		}
		std::cout << "Detected " << group.getNumObjs() << " variables in file.\n";*/

		H5::Group group_meta(file.openGroup("universe"));
		/*for (int grps = 0; grps < group_meta.getNumAttrs(); grps++)
		{
			H5::Attribute attr = group_meta.openAttribute(grps);
			std::cout << "Universe: " << attr.getName() << " : ";
			std::string name = attr.getName();
			double val = 0.0;
			H5::DataType type = attr.getDataType();
			attr.read(type, &val);
			std::cout << val << std::endl;
		}
		std::cout << "Detected " << group_meta.getNumAttrs() << " universe attributes.\n";*/

		int fields = group.getNumObjs();

		H5::DataSet dataset(group.openDataSet(paramName));
		H5::DataSpace dataspace(dataset.getSpace());
		H5::DataSpace memspace(dataset.getSpace()); //This would define rank and local rank extent
		hsize_t tdims[3];
		dataspace.getSimpleExtentDims(tdims);
		//std::cout << "Data dimensions: " << dims[0] << " " << dims[1] << " " << dims[2] << "\n";
		numElements = tdims[0] * tdims[1] * tdims[2];
		dims[0] = tdims[0];
		dims[1] = tdims[1];
		dims[2] = tdims[2];
			
		totalNumberOfElements = numElements; // Temporary 

		dataType = "float";
		// Set-up data stream
		allocateMem(dataType, numElements, 0);

		// Data buffer stream, data_type, memory_space, file_space
		// H%::PredType::NATIVE_DOUBLE
		dataset.read(data, H5::PredType::NATIVE_FLOAT, memspace, dataspace);

		dataset.close();
		file.close();
	}
	catch (H5::FileIException error)
	{
		//error.printError();
		return -1;
	}
	// catch failure caused by the DataSet operations
	catch (H5::DataSetIException error)
	{
		//error.printError();
		return -1;
	}
	// catch failure caused by the DataSpace operations
	catch (H5::DataSpaceIException error)
	{
		//error.printError();
		return -1;
	}
	// catch failure caused by the DataSpace operations
	catch (H5::DataTypeIException error)
	{
		//error.printError();
		return -1;
	}


	clock.stop();
	log << "Loading data took " << clock.getDuration() << " s" << std::endl;

	return 1; // All good
}

inline int NYXDataLoader::saveCompData(std::string paramName, void * cData)
{
	compFullData.insert({ paramName, cData });
	return 1;
}

inline int NYXDataLoader::writeData(std::string _filename)
{
	return 1;
}

#endif
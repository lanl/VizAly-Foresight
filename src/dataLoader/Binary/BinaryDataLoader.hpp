/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Jesus Pulido
================================================================================*/

#ifndef _BINARY_LOADER_H_
#define _BINARY_LOADER_H_

#include <sstream>
#include <string>
#include "dataLoaderInterface.hpp"

// Helper Functions
#include "json.hpp"
#include "timer.hpp"
#include <unordered_map>

class BinaryDataLoader: public DataLoaderInterface
{
	int numRanks;
	int myRank;
	size_t headerSize;
	std::string dataTarget;

	std::unordered_map<std::string, void *> compFullData;

  public:
	BinaryDataLoader();
	~BinaryDataLoader();
	int allocateMem(std::string dataType, size_t numElements, int offset);
	int deAllocateMem();
	size_t findByteAddress(size_t x, size_t y, size_t z, size_t field);

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int saveCompData(std::string paramName, void * cData);
	int writeData(std::string _filename);
    int saveInputFileParameters() { return 1; };
	int close() { return deAllocateMem(); }
	void setParam(std::string paramName, std::string type, std::string value){};
};


inline BinaryDataLoader::BinaryDataLoader()
{
	myRank = 0;
	numRanks = 0;
	loader = "Binary";
	saveData = false;
}

inline BinaryDataLoader::~BinaryDataLoader()
{
	deAllocateMem();
}


inline void BinaryDataLoader::init(std::string _filename, MPI_Comm _comm)
{
	filename = _filename;
	comm = _comm;
	headerSize = 0;
	saveData = false;

	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);

	std::unordered_map<std::string, std::string>::const_iterator got = loaderParams.find("dims");

	// parse and set dims
	if (got != loaderParams.end())
	{
		int cnt = 0;
		size_t pos = 0;
		std::string token; std::string s = loaderParams["dims"];

		// delete "[" and "]"
		s.erase(0, 1); s.erase(s.size()-1, s.size());

		// extract dims using delimiter
		while ((pos = s.find(",")) != std::string::npos) {
			token = s.substr(0, pos);
			sizePerDim[cnt] = strConvert::to_uint64(token);
			s.erase(0, pos + 1);
			cnt++;
		}
		sizePerDim[cnt] = strConvert::to_uint64(token);
	}

	got = loaderParams.find("header");
	if (got != loaderParams.end())
	{
		headerSize = strConvert::to_uint64(got->second);
	}

	got = loaderParams.find("type");
	if (got != loaderParams.end())
	{
		std::string s = got->second;
		// delete " and "
		s.erase(0, 1); s.erase(s.size() - 1, s.size());

		dataType = s;
	}

	got = loaderParams.find("target");
	if (got != loaderParams.end())
	{
		std::string s = got->second;
		// delete " and "
		s.erase(0, 1); s.erase(s.size() - 1, s.size());

		dataTarget = s;
	}

	std::cout << "sizePerDim: " << sizePerDim[0] << "," << sizePerDim[1] << "," << sizePerDim[2] << ", header: " << headerSize << ", type: " << dataType << std::endl;
}

inline int BinaryDataLoader::allocateMem(std::string dataType, size_t numElements, int offset)
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

inline int BinaryDataLoader::deAllocateMem()
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

inline int BinaryDataLoader::loadData(std::string paramName)
{
	Timer clock;
	log.str("");

	clock.start();
	param = paramName;

	int field = 0;
	field = strConvert::to_int(paramName);

	// Set Dims
	size_t dimx=sizePerDim[0];
	size_t dimy=sizePerDim[1];
	size_t dimz=sizePerDim[2];

	totalNumberOfElements = dimx * dimy * dimz;
	numElements = totalNumberOfElements;

	std::ifstream f;

	// Try opening the file
	f.open(filename, std::ifstream::binary);
	if (f.fail()) {
		std::cout << "Error opening file" << std::endl;
		return -1;
	}

	f.seekg(0, std::ios::end);
	double file_len = f.tellg();
	f.seekg(0, std::ios::beg);

	size_t field_len = dimx*dimy*dimz*elemSize;

	size_t xl = 0; size_t yl = 0; size_t zl = 0;
	size_t xu = dimx - 1; size_t yu = dimy - 1; size_t zu = dimz - 1;

	allocateMem(dataTarget, numElements, 0);

	float * fdata; double * ddata;
	
	size_t cnt = 0;
	if (dataTarget == "float")
	{
		fdata = static_cast<float*>(data);
		for (size_t z = zl; z <= zu; z++)
		{
			for (size_t y = yl; y <= yu; y++)
			{
				std::streampos seekVal = findByteAddress(xl, y, z, field);
				f.seekg(seekVal);

				if (dataType == "float")
				{
					float value;
					for (size_t x = xl; x <= xu; x++)
					{
						f.read(reinterpret_cast<char*>(&value), sizeof(float));
						fdata[cnt] = value; cnt++;
					}
				}
				else if (dataType == "double")
				{
					double value;
					for (size_t x = xl; x <= xu; x++)
					{
						f.read(reinterpret_cast<char*>(&value), sizeof(double));
						fdata[cnt] = value; cnt++;
					}
				}
				else {}
			}
		}
	}
	else if (dataTarget == "double")
	{
		ddata = static_cast<double*>(data);
		for (size_t z = zl; z <= zu; z++)
		{
			for (size_t y = yl; y <= yu; y++)
			{
				std::streampos seekVal = findByteAddress(xl, y, z, field);
				f.seekg(seekVal);

				if (dataType == "float")
				{
					float value;
					for (size_t x = xl; x <= xu; x++)
					{
						f.read(reinterpret_cast<char*>(&value), sizeof(float));
						ddata[cnt] = value; cnt++;
					}
				}
				else if (dataType == "double")
				{
					double value;
					for (size_t x = xl; x <= xu; x++)
					{
						f.read(reinterpret_cast<char*>(&value), sizeof(double));
						ddata[cnt] = value; cnt++;
					}
				}
				else {}
			}
		}
	}
	else
	{
		std::cout << "Error: Unsupported binary dataType!" << std::endl; return -1;
	}


	f.close();

	numElements = cnt;

	clock.stop();
	log << "Loading data took " << clock.getDuration() << " s" << std::endl;

	return 1; // All good
}

inline int BinaryDataLoader::saveCompData(std::string paramName, void * cData)
{
	compFullData.insert({ paramName, cData });
	return 1;
}

inline int BinaryDataLoader::writeData(std::string _filename)
{
	return 1;
}

// Given x,y,z coordinates and current scalar field, this will find the
// byte location of the data within a binary file
size_t BinaryDataLoader::findByteAddress(size_t x, size_t y, size_t z, size_t field)
{
	// add the header
	size_t byteLocation = headerSize;

	// first look at the scalar_field
	byteLocation += field*sizePerDim[0]*sizePerDim[1]*sizePerDim[2]*elemSize;

	// Look at z
	byteLocation += z*sizePerDim[0]*sizePerDim[1]*elemSize;

	// look at y
	byteLocation += y*sizePerDim[0]*elemSize;

	// look at x
	byteLocation += x*elemSize;

	return byteLocation;
}

#endif
/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/

#ifndef _DATA_LOADER_INTERFACE_H_
#define _DATA_LOADER_INTERFACE_H_

#include <string>
#include <mpi.h>
#include <sstream>
#include <unordered_map>
#include <list>

#include "HACC/gioData.hpp"
#include "log.hpp"


class DataLoaderInterface 
{
  protected:
	std::string loader;
	std::string filename;
	std::string dataType;
	std::string param;

	bool saveData;
	int origNumDims;
	int timestep;


	size_t origDims[5]{ 0,0,0,0,0 };	// Global dataset size
	size_t sizePerDim[5]{ 0,0,0,0,0 };	// For compression, size of this mpi rank
	size_t rankOffset[3]{0, 0, 0};		// Rank offset, for parallel operations
	size_t elemSize;					// sizeof() in bytes of that parameter
	size_t totalNumberOfElements;		// total number of points/cells for input file
	size_t numElements;					// number of points/cells for this mpi rank

	MPI_Comm comm;						// global mpi handler

  public:   // TO_CHANGE, *data should not be public
	void *data;							// raw data pointer (do not touch outside!!)
	
	std::unordered_map<std::string, std::string> loaderParams;  // data-specific parameters, used by binary reader
	std::vector<GioData> inOutData;				// passthrough data, used by genericio reader
	std::vector<std::string> inOutDataName;		// passthrough data

  public:
	virtual void init(std::string _filename, MPI_Comm _comm) = 0;
	virtual int loadData(std::string paramName) = 0;					// Read op
	virtual int saveCompData(std::string paramName, void * cData) = 0;	// Stores the data-pointer for future writing
	virtual int writeData(std::string _filename) = 0;					// Write op
	virtual int saveInputFileParameters() = 0;							// Reads file header only and stores data properties
	virtual int close() = 0;
	virtual void setParam(std::string paramName, std::string type, std::string value) = 0;
	virtual bool loadUncompressedFields(nlohmann::json const& jsonInput) = 0; // Data passthrough, used by hdf5 reader

	void setTimestep(int ts){ timestep = ts; }
	size_t getNumElements() { return numElements; }
	size_t * getSizePerDim() { return sizePerDim; }
	size_t getTypeSize() { return elemSize; }
	std::string getType() { return dataType; }
	std::string getParam() { return param; }


	void setSave(bool state) { saveData = state; }	// if true, write out decomp data


	std::string getDataInfo()
	{
		std::stringstream dataInfo;
		dataInfo << "\nLoader type: " << loader << std::endl;
		dataInfo << "Filename: " << filename << std::endl;
		dataInfo << "Total number of elements: " << totalNumberOfElements << std::endl;
		dataInfo << "Param: " << param << std::endl;
		dataInfo << "dataType: " << dataType << std::endl;
		dataInfo << "numElements: " << numElements << std::endl;
		dataInfo << "sizePerDim: " << sizePerDim[0] << " " << sizePerDim[1] 
           << " " << sizePerDim[2] << " " << sizePerDim[3] << " " 
           << sizePerDim[4] << std::endl;

		return dataInfo.str();
	}
};



class Partition
{
  public:
	int min_x = 0;
	int min_y = 0;
	int min_z = 0;
	int max_x = 0;
	int max_y = 0;
	int max_z = 0;

  public:
	Partition() = default;
	Partition(
		int in_min_x, int in_min_y, int in_min_z,
		int in_max_x, int in_max_y, int in_max_z
	) : min_x(in_min_x),
		min_y(in_min_y),
		min_z(in_min_z),
		max_x(in_max_x),
		max_y(in_max_y),
		max_z(in_max_z)
	{
	}

	void print()
	{
		std::cout << "Dims: " << min_x << ", " << min_y << ", " << min_z << " - "
			<< max_x << ", " << max_y << ", " << max_z << std::endl;
	}
};

inline Partition getPartition(int myRank, int numRanks, int extentsX, int extentsY, int extentsZ)
{
	std::list<Partition> partitions;
	partitions.push_back(Partition{ 0, 0, 0, extentsX, extentsY, extentsZ });


	Partition first_half, second_half;

	int axis = 0;
	while (partitions.size() < numRanks)
	{
		int numCurrentPartitions = partitions.size();
		for (int i = 0; i < numCurrentPartitions; i++)
		{
			Partition parent = partitions.front();
			partitions.pop_front();


			if (axis == 0)
			{
				if (parent.max_x - parent.min_x <= 1)
				{
					partitions.push_back(parent);
					continue;
				}

				first_half.min_x = parent.min_x;
				first_half.max_x = (parent.min_x + parent.max_x) / 2;
				second_half.min_x = first_half.max_x;
				second_half.max_x = parent.max_x;

				first_half.min_y = second_half.min_y = parent.min_y;
				first_half.max_y = second_half.max_y = parent.max_y;
				first_half.min_z = second_half.min_z = parent.min_z;
				first_half.max_z = second_half.max_z = parent.max_z;
			}

			if (axis == 1)
			{
				if (parent.max_y - parent.min_y <= 1)
				{
					partitions.push_back(parent);
					continue;
				}

				first_half.min_y = parent.min_y;
				first_half.max_y = (parent.min_y + parent.max_y) / 2;
				second_half.min_y = first_half.max_y;
				second_half.max_y = parent.max_y;

				first_half.min_x = second_half.min_x = parent.min_x;
				first_half.max_x = second_half.max_x = parent.max_x;
				first_half.min_z = second_half.min_z = parent.min_z;
				first_half.max_z = second_half.max_z = parent.max_z;
			}

			if (axis == 2)
			{
				if (parent.max_z - parent.min_z <= 1)
				{
					partitions.push_back(parent);
					continue;
				}

				first_half.min_z = parent.min_z;
				first_half.max_z = (parent.min_z + parent.max_z) / 2;
				second_half.min_z = first_half.max_z;
				second_half.max_z = parent.max_z;

				first_half.min_x = second_half.min_x = parent.min_x;
				first_half.max_x = second_half.max_x = parent.max_x;
				first_half.min_y = second_half.min_y = parent.min_y;
				first_half.max_y = second_half.max_y = parent.max_y;
			}

			partitions.push_back(first_half);
			partitions.push_back(second_half);

			if (partitions.size() >= numRanks)
				break;
		}

		axis++;
		if (axis == 3)
			axis = 0;
	}


	auto it = partitions.begin();
	std::advance(it, myRank);

	return *it;
}


inline void getMPIDivisions(int numRanks, int numDims, int divisions[3])
{
	divisions[0] = divisions[1] = divisions[2] = 1;
	
	int axis = 0;
	while (divisions[0]*divisions[1]*divisions[2] < numRanks)
	{
		if (axis == 0)
			divisions[0] = divisions[0]*2;
		else if (axis == 1)
			divisions[1] = divisions[1]*2;
		else if (numDims > 2)
			if (axis == 2)
				divisions[2] = divisions[2]*2;

		axis = axis + 1;
		if (axis == numDims)
			axis = 0;
	}
}

#endif

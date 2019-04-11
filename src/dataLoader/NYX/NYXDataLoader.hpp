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

#include <list>
#include <iterator>
#include <algorithm>


struct partition
{ 
	int min_x; 
	int min_y; 
	int min_z;

	int max_x; 
	int max_y; 
	int max_z;

	partition(){};
	partition(int _min_x, int _min_y, int _min_z, int _max_x, int _max_y, int _max_z){ setPatition(_min_x,_min_y,_min_z, _max_x,_max_y,_max_z); }
	void setPatition(int _min_x, int _min_y, int _min_z, int _max_x, int _max_y, int _max_z)
	{ 
		min_x = _min_x; 	max_x = _max_x; 
		min_y = _min_y; 	max_y = _max_y; 
		min_z = _min_z; 	max_z = _max_z; 
	}

	void print()
	{
		std::cout << min_x << ", " << min_y << ", " << min_z << " - "
				  << max_x << ", " << max_y << ", " << max_z << std::endl;
	}
};


partition getPartition(int myRank, int numRanks, int extentsX, int extentsY, int extentsZ)
{
	float extents[3];
	extents[0] = extentsX;
	extents[1] = extentsY;
	extents[2] = extentsZ;

	std::list<partition> partitions;
	partition temp(0, 0, 0, extentsX, extentsY, extentsZ);
	partitions.push_back(temp);

	partition firstHalf, secondHalf;

	int axis = 0;
	while (partitions.size() < numRanks)
	{
		int numCurrentPartitions = partitions.size();
		for (int i=0; i<numCurrentPartitions; i++)
		{
			partition parent = partitions.front();
			partitions.pop_front();

			// if (myRank == 0)
			// {
			// 	std::cout << "axis: " << axis << std::endl; 
			// 	std::cout << "Parent: "; 
			// 	parent.print();
			// }


			if (axis == 0)
			{
				firstHalf.min_x = parent.min_x;
				firstHalf.max_x = (parent.min_x + parent.max_x)/2;

				secondHalf.min_x = firstHalf.max_x;
				secondHalf.max_x = parent.max_x;

				firstHalf.min_y = secondHalf.min_y = parent.min_y;
				firstHalf.max_y = secondHalf.max_y = parent.max_y;

				firstHalf.min_z = secondHalf.min_z = parent.min_z;
				firstHalf.max_z = secondHalf.max_z = parent.max_z;
			}

			if (axis == 1)
			{
				firstHalf.min_y = parent.min_y;
				firstHalf.max_y = (parent.min_y + parent.max_y)/2;

				secondHalf.min_y = firstHalf.max_y;
				secondHalf.max_y = parent.max_y;

				firstHalf.min_x = secondHalf.min_x = parent.min_x;
				firstHalf.max_x = secondHalf.max_x = parent.max_x;

				firstHalf.min_z = secondHalf.min_z = parent.min_z;
				firstHalf.max_z = secondHalf.max_z = parent.max_z;
			}

			if (axis == 2)
			{
				firstHalf.min_z = parent.min_z;
				firstHalf.max_z = (parent.min_z + parent.max_z)/2;

				secondHalf.min_z = firstHalf.max_z;
				secondHalf.max_z = parent.max_z;

				firstHalf.min_x = secondHalf.min_x = parent.min_x;
				firstHalf.max_x = secondHalf.max_x = parent.max_x;

				firstHalf.min_y = secondHalf.min_y = parent.min_y;
				firstHalf.max_y = secondHalf.max_y = parent.max_y;
			}

			partitions.push_back(firstHalf);
			partitions.push_back(secondHalf);


			// if (myRank == 0)
			// {
			// 	std::cout << "child 1: ";	firstHalf.print();
			// 	std::cout << "child 2: ";	secondHalf.print();
			// 	std::cout << "\n";
			// }

			if (partitions.size() >= numRanks)
				break;
		}


		axis++;
			if (axis == 3)
				axis = 0;



		// if (myRank == 0)
		// {
		// 	for (auto it=partitions.begin(); it!=partitions.end() ; it++)
		// 	{
		// 		(*it).print();
		// 	}

		// 	std::cout << "\n\n";
		// }
	}


	auto it = partitions.begin();
	std::advance(it, myRank);

	return *it;

}  	


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
		H5::Group group_meta(file.openGroup("universe"));

		int fields = group.getNumObjs();

		H5::DataSet dataset(group.openDataSet(paramName));
		H5::DataSpace dataspace(dataset.getSpace());
		H5::DataSpace memspace(dataset.getSpace()); //This would define rank and local rank extent
		hsize_t tdims[3];
		dataspace.getSimpleExtentDims(tdims);
		totalNumberOfElements = tdims[0] * tdims[1] * tdims[2];

		log << "Param: " << paramName << std::endl;
		log << "Data dimensions: " << tdims[0] << " " << tdims[1] << " " << tdims[2] <<  " | totalNumberOfElements " << totalNumberOfElements << "\n";
		

		// Read only a subset of the file
		hsize_t count[3];              
    	hsize_t offset[3];            


    	partition current = getPartition(myRank, numRanks, tdims[0], tdims[1], tdims[2]);
    	count[0] = current.max_x - current.min_x;
    	count[1] = current.max_y - current.min_y;
    	count[2] = current.max_z - current.min_z;

    	offset[0] = current.min_x;
    	offset[1] = current.min_y;
    	offset[2] = current.min_z;

		numElements = count[0] * count[1] * count[2];


        log << myRank << " ~ Count : " << count[0] << " " << count[1] << " " << count[2]  << " | " << numElements << "\n";
		log << myRank << " ~ Offset: " << offset[0] << " " << offset[1] << " " << offset[2] << "\n";

        dataspace.selectHyperslab( H5S_SELECT_SET, count, offset );
        H5::DataSpace memspaceRead( 3, count, NULL );
        hsize_t zero[3] = {0,0,0};
        memspaceRead.selectHyperslab( H5S_SELECT_SET, count, zero );
       

        // Set-up data stream
		dataType = "float";
		allocateMem(dataType, numElements, 0);


		dataset.read(data, H5::PredType::NATIVE_FLOAT, memspaceRead, dataspace);


		for (int z=0; z<count[2]; z++)
		{
			for (int y=0; y<count[1]; y++)
			{
				for (int x=0; x<count[0]; x++)
				{
					log << ((float *)data)[z*count[0]*count[1] + y*count[0]+ x] << " ";
				}
				log << std::endl;
			}

			log << std::endl << std::endl;
		}
     

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
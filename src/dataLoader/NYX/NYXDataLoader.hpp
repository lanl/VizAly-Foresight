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
//#include "H5Cpp.h"
#include "hdf5.h"

#include <list>
#include <iterator>
#include <algorithm>


struct uncompressedData
{
	std::string paramName;
	std::string dataType;
	size_t numElements;
	size_t dims[3];
	void *data;

	uncompressedData(){};

	uncompressedData(std::string _paramName, std::string _dataType, size_t _numElements, size_t _dims[3])
	{
		paramName = _paramName;
		dataType = _dataType;
		numElements = _numElements;

		for (int i=0; i<3; i++)
			dims[i] = _dims[i];
	}

	inline int allocateMem()
	{
	    if (dataType == "float")
	        data = new float[numElements];
	    else if (dataType == "double")
	        data = new double[numElements];
	    else if (dataType == "int8_t")
	        data = new int8_t[numElements];
	    else if (dataType == "int16_t")
	        data = new int16_t[numElements];
	    else if (dataType == "int32_t")
	        data = new int32_t[numElements];
	    else if (dataType == "int64_t")
	        data = new int64_t[numElements];
	    else if (dataType == "uint8_t")
	        data = new uint8_t[numElements];
	    else if (dataType == "uint16_t")
	        data = new uint16_t[numElements];
	    else if (dataType == "uint32_t")
	        data = new uint32_t[numElements];
	    else if (dataType == "uint64_t")
	        data = new uint64_t[numElements];
	    else
	        return 0;

	    return 1;
	}


	inline int deAllocateMem()
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
};


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

	std::vector<uncompressedData> toWriteData;

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

	inOutData.clear();
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
	//try {
	{
		log << "filename: " << filename << " param name: " << paramName << std::endl;

		hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		hid_t group = H5Gopen(file, "native_fields", H5P_DEFAULT);
		hid_t group_meta = H5Gopen(file, "universe", H5P_DEFAULT);

		hsize_t fields;
		herr_t err = H5Gget_num_objs(group, &fields);


		paramName = "/native_fields/" + paramName;
		log << "paramName: " <<  paramName << std::endl;
		hid_t dataset = H5Dopen(file, paramName.c_str(), H5P_DEFAULT);
		hid_t dataspace = H5Dget_space(dataset);
		hid_t memspace = H5Dget_space(dataset);


		hsize_t tdims[3];
		int ndims = H5Sget_simple_extent_dims(dataspace, tdims, NULL);
		origDims[0] = tdims[0];
		origDims[1] = tdims[1];
		origDims[2] = tdims[2];

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

    	sizePerDim[0] = count[0];
		sizePerDim[1] = count[1];
		sizePerDim[2] = count[2];

    	offset[0] = current.min_x;
    	offset[1] = current.min_y;
    	offset[2] = current.min_z;

		numElements = count[0] * count[1] * count[2];

		rankOffset[0] = offset[0];
		rankOffset[1] = offset[1];
		rankOffset[2] = offset[2];


        log << myRank << " ~ Count : " << count[0]  << " " << count[1]  << " " << count[2]  << " | " << numElements << "\n";
		log << myRank << " ~ Offset: " << offset[0] << " " << offset[1] << " " << offset[2] << "\n";


		herr_t status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
        hid_t memspaceRead = H5Screate_simple(3, count, NULL);

        hsize_t zero[3] = {0,0,0};
        status = H5Sselect_hyperslab(memspaceRead, H5S_SELECT_SET, zero, NULL, count, NULL);
       

        // Set-up data stream
		dataType = "float";
		allocateMem(dataType, numElements, 0);
		log << myRank << " ~ numElements: " <<numElements << "\n";

		status = H5Dread(dataset, H5T_NATIVE_FLOAT, memspaceRead, dataspace, H5P_DEFAULT, data);


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


		log << myRank << " ~ done!!! " << "\n";
     

		H5Dclose(dataset);
		H5Fclose(file);
	}
	
	clock.stop();
	log << "Loading data took " << clock.getDuration() << " s" << std::endl;

	return 1; // All good
}



inline int NYXDataLoader::saveCompData(std::string paramName, void * cData)
{
	size_t _dims[3];
	_dims[0] = sizePerDim[0];
	_dims[1] = sizePerDim[1];
	_dims[2] = sizePerDim[2];

	uncompressedData temp(paramName, "float", numElements, _dims);
	temp.allocateMem();
	std::memcpy(temp.data, cData, numElements*sizeof(float));

	compFullData.insert({ paramName, temp.data });
	toWriteData.push_back(temp);

	return 1;
}

inline int NYXDataLoader::writeData(std::string _filename)
{
	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);

    //std::cout << "_filename: " << _filename << std::endl;


	hid_t file_id = H5Fcreate(_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);

	hsize_t fileDims[3];
	fileDims[0] = origDims[0];
	fileDims[1] = origDims[1];
	fileDims[2] = origDims[2];
	hid_t filespace = H5Screate_simple(3, fileDims, NULL);

	hid_t group = H5Gcreate2(file_id, "/native_fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	int numDatasets = compFullData.size();
	hid_t *dset_id = new hid_t[numDatasets];

	int datasetCount = 0;
	for (auto& scalar: compFullData)
	{
		std::string fieldname = "/native_fields/" + scalar.first;
		//std::cout << myRank << " ~ inOutData[i].name: " << scalar.first << std::endl;

		dset_id[datasetCount] = H5Dcreate(file_id, fieldname.c_str(), H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    datasetCount++;
	}

    hsize_t	 count[3];	          
    count[0] = sizePerDim[0];
    count[1] = sizePerDim[1];
    count[2] = sizePerDim[2];

    hsize_t	offset[3];
    offset[0] = rankOffset[0];
    offset[1] = rankOffset[1];
    offset[2] = rankOffset[2];
    hid_t memspace = H5Screate_simple(3, count, NULL);


    filespace = H5Dget_space(dset_id[0]);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    datasetCount = 0;
    for (auto& scalar: compFullData)
	{
		log << "outputting " << scalar.first << std::endl;
		for (int z=0; z<count[2]; z++)
		{
			for (int y=0; y<count[1]; y++)
			{
				for (int x=0; x<count[0]; x++)
				{
					log << ((float *)scalar.second)[z*count[0]*count[1] + y*count[0]+ x] << " ";
				}
				log << std::endl;
			}

			log << std::endl << std::endl;
		}


	    herr_t status = H5Dwrite(dset_id[datasetCount], H5T_NATIVE_FLOAT, memspace, filespace, plist_id, (float*)scalar.second);
	    datasetCount++;
	}

	for (int i=0; i<numDatasets; i++)
    	H5Dclose(dset_id[i]);

	//H5Dclose(filespace);
	//H5Sclose(memspace);


	herr_t status = H5Gclose(group);


    H5Pclose(plist_id);
    H5Fclose(file_id);


    for (int i=0; i<toWriteData.size(); i++)
    	toWriteData[i].deAllocateMem();
}

#endif
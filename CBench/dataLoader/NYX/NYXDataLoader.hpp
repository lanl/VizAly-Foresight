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
#include "hdf5.h"

#include <list>
#include <iterator>
#include <algorithm>

class UncompressedData{

  public:
	std::string paramName = "";
	std::string dataType = "";
	size_t numElements = 0;
	size_t dims[3] = {0, 0, 0};
	void *data = nullptr;

  public:
	UncompressedData() = default;
	UncompressedData(UncompressedData const&) = default;
	UncompressedData(UncompressedData&&) = default;
	UncompressedData(
		std::string in_param_name,
		std::string in_data_type,
		size_t in_numElements,
		size_t in_dims[3]
	) : paramName(in_param_name),
		dataType (in_data_type),
		numElements (in_numElements),
		dims { in_dims[0], in_dims[1], in_dims[2] }
	{}

	inline bool allocateMem()
	{
		return Memory::allocate(data, dataType, numElements);
	}

	inline bool deAllocateMem()
	{
		return Memory::release(data, dataType);
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
	{}

	void print()
	{
		std::cout << min_x << ", " << min_y << ", " << min_z << " - "
				  << max_x << ", " << max_y << ", " << max_z << std::endl;
	}
};

Partition getPartition(int myRank, int numRanks, int extentsX, int extentsY, int extentsZ)
{

	std::list<Partition> partitions;
	partitions.push_back(Partition{0, 0, 0, extentsX, extentsY, extentsZ});

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


class NYXDataLoader : public DataLoaderInterface
{
	int myRank = 0;
	int numRanks = 0;
	std::string defaultGroup = "native_fields";

	std::vector<UncompressedData> toWriteData;
	// group data: key=group, value:array of field data
	std::unordered_map<std::string, std::vector<UncompressedData>> groupsData;
	std::unordered_map<std::string, std::vector<UncompressedData>> groupsAttribs;

  public:
	NYXDataLoader() { loader = "NYX"; saveData = false; }
	~NYXDataLoader() { deAllocateMem(); }

	bool allocateMem(std::string dataType, size_t numElements, int offset);
	bool deAllocateMem();

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int loadDataOld(std::string paramName);
	int saveCompData(std::string paramName, void *cData);
	int writeData(std::string _filename);
	int writeDataOld(std::string _filename);
	int saveInputFileParameters() { return 1; }
	int close() { return deAllocateMem(); }
	void setParam(std::string paramName, std::string type, std::string value);
	bool loadUncompressedFields(nlohmann::json const& jsonInput);

  private:
	hid_t getNativeDataType(std::string const& dataType) const;

	void writeGroupAttrib(
		hid_t& in_file_id,
		std::string const& in_group_name
	);

	void writeGroupData(
		hid_t& in_file_id, hid_t& in_filespace, hid_t& in_memspace,
		std::string const& in_group_name,
		hsize_t const count[3], hsize_t const offset[3]
	);

	herr_t loadGroupDataSet(
		hid_t const& in_file,
		std::string const& in_group_name,
		std::string const& in_field_name
	);

	herr_t loadAttribute(
		hid_t const& in_file, std::string const& in_group_name,
		std::string const& in_attrib_name, int num_elems = 1
	);
};

inline void NYXDataLoader::init(std::string _filename, MPI_Comm _comm)
{
	filename = _filename;
	comm = _comm;
	saveData = false;

	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);
}

inline void NYXDataLoader::setParam(std::string paramName, std::string type, std::string value)
{
	if (paramName == "group")
		defaultGroup = value;
}

// TODO: refactor it in a util class
inline bool NYXDataLoader::allocateMem(std::string dataType, size_t numElements, int offset)
{
	elemSize = Memory::sizeOf[dataType];
	return Memory::allocate(data, dataType, numElements, offset);
}


inline bool NYXDataLoader::deAllocateMem()
{
	return Memory::release(data, dataType);
}

inline int NYXDataLoader::loadData(std::string paramName)
{

	Timer clock;
	clock.start();

	log.str("");
	log << "filename: " << filename << " param name: " << paramName << std::endl;
	hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	herr_t status = loadGroupDataSet(file, defaultGroup, paramName);

	H5Fclose(file);
	clock.stop();
	log << "Loading data took " << clock.getDuration() << " s" << std::endl;

	int const exit_code = (status < 0 ? 0 : 1);
	return exit_code;
}


// generic group data field loader
// assume that all datasets have the same dimension
inline herr_t NYXDataLoader::loadGroupDataSet(
	hid_t const& in_file,
	std::string const& in_group_name,
	std::string const& in_field_name
)
{

	hsize_t fields;
	hsize_t tdims[3];
	hsize_t count[3];
	hsize_t offset[3];
	hsize_t zero[3] = {0, 0, 0};

	std::string param_name = "/" + in_group_name + "/" + in_field_name;
	log << "paramName: " << param_name << std::endl;

	hid_t group     = H5Gopen(in_file, in_group_name.c_str(), H5P_DEFAULT);
	herr_t err      = H5Gget_num_objs(group, &fields);
	hid_t dataset   = H5Dopen(in_file, param_name.c_str(), H5P_DEFAULT);
	hid_t dataspace = H5Dget_space(dataset);

	int ndims = H5Sget_simple_extent_dims(dataspace, tdims, NULL);
//  log << "origDims:("<< origDims[0] <<", "<< origDims[1] <<", "<< origDims[2] <<")"<< std::endl;
//  log << "tdims:("<< (unsigned) tdims[0] <<", "<< (unsigned) tdims[1] <<", "<< (unsigned) tdims[2] <<")"<< std::endl;

	origDims[0] = tdims[0];
	origDims[1] = tdims[1];
	origDims[2] = tdims[2];

	totalNumberOfElements = tdims[0] * tdims[1] * tdims[2];

	log << "Param: " << param_name << std::endl;
	log << "Data dimensions: " << tdims[0] << " " << tdims[1] << " " << tdims[2]
		<< " | totalNumberOfElements " << totalNumberOfElements << std::endl;

	// Read only a subset of the file
	Partition current = getPartition(myRank, numRanks, tdims[0], tdims[1], tdims[2]);
	count[0] = current.max_x - current.min_x;
	count[1] = current.max_y - current.min_y;
	count[2] = current.max_z - current.min_z;

	sizePerDim[0] = count[0];
	sizePerDim[1] = count[1];
	sizePerDim[2] = count[2];

	numElements = count[0] * count[1] * count[2];

	offset[0] = current.min_x;
	offset[1] = current.min_y;
	offset[2] = current.min_z;

	rankOffset[0] = offset[0];
	rankOffset[1] = offset[1];
	rankOffset[2] = offset[2];

	log << myRank << " ~ Count : " << count[0] << " " << count[1] << " " << count[2] << " | " << numElements << "\n";
	log << myRank << " ~ Offset: " << offset[0] << " " << offset[1] << " " << offset[2] << "\n";
	log << myRank << " ~ numElements: " << numElements << "\n";

	// select data section
	hid_t memspace = H5Screate_simple(3, count, NULL);
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
	H5Sselect_hyperslab(memspace, H5S_SELECT_SET, zero, NULL, count, NULL);

	// Set-up data stream and select only float values
	bool const is_default = (in_group_name == defaultGroup);
	// two cases: default group and custom group
	if (is_default)
	{
		dataType = "float";
		allocateMem(dataType, numElements, 0);
	}
	else
	{
		auto& group_data = groupsData[in_group_name];
		group_data.emplace_back(in_field_name, "float", numElements, sizePerDim);
		auto& field_dataset = group_data.back();
		field_dataset.allocateMem();
	}

	// set data storage pointer and copy dataset
	auto storage = (is_default ? data : groupsData[in_group_name].back().data);
	hid_t const native_datatype = getNativeDataType(dataType);
	herr_t status = H5Dread(dataset, native_datatype, memspace, dataspace, H5P_DEFAULT, storage);

	H5Dclose(dataset);
	return status;
}

// (!) attributes may be scalars or arrays of any datatype
inline herr_t NYXDataLoader::loadAttribute(
	hid_t const& in_file, std::string const& in_group_name,
	std::string const& in_attrib_name, int num_elems
)
{

	assert(num_elems > 0);
	std::string param_name = "/" + in_group_name + "/" + in_attrib_name;
	log << "attrib_name: " << param_name << std::endl;

	hid_t group_id  = H5Gopen(in_file, in_group_name.c_str(), H5P_DEFAULT);
	hid_t attrib_id = H5Aopen(group_id, in_attrib_name.c_str(), H5P_DEFAULT);
	hid_t datatype  = H5Aget_type(attrib_id);

	std::string type = "";
	size_t size_per_dim[] = {1, 1, 1};  // arbitrary

	switch (H5Tget_class(datatype))
	{
	case H5T_INTEGER: type = "int"; break;
	case H5T_FLOAT:   type = "double"; break;
	case H5T_STRING:  type = "string"; break;
	default: throw std::runtime_error("Bad group attribute type");
	}

	// set up and uncompressed data object
	auto& group_data = groupsAttribs[in_group_name];
	group_data.emplace_back(in_attrib_name, type, num_elems, size_per_dim);
	auto& current_attrib = group_data.back();
	assert(current_attrib.allocateMem());
	assert(current_attrib.data != nullptr);
	// read and store the group attribute
	herr_t status = H5Aread(attrib_id, datatype, current_attrib.data);

	H5Tclose(datatype);
	H5Aclose(attrib_id);
	H5Gclose(group_id);
	return status;
}

inline int NYXDataLoader::saveCompData(std::string paramName, void *cData)
{
	size_t _dims[3];
	_dims[0] = sizePerDim[0];
	_dims[1] = sizePerDim[1];
	_dims[2] = sizePerDim[2];

	UncompressedData temp(paramName, "float", numElements, _dims);
	temp.allocateMem();
	std::memcpy(temp.data, cData, numElements * sizeof(float));

	toWriteData.push_back(temp);

	return 1;
}

// TODO: replace it by a lookup table
inline hid_t NYXDataLoader::getNativeDataType(std::string const& dataType) const
{

	if (dataType == "int")
	{
		return H5T_NATIVE_INT;
	}
	else if (dataType == "float")
	{
		return H5T_NATIVE_FLOAT;
	}
	else if (dataType == "double")
	{
		return H5T_NATIVE_DOUBLE;
	}
	else if (dataType == "int8_t")
	{
		return H5T_NATIVE_INT8;
	}
	else if (dataType == "int16_t")
	{
		return H5T_NATIVE_INT16;
	}
	else if (dataType == "int32_t")
	{
		return H5T_NATIVE_INT32;
	}
	else if (dataType == "int64_t")
	{
		return H5T_NATIVE_INT64;
	}
	else if (dataType == "uint8_t")
	{
		return H5T_NATIVE_UINT8;
	}
	else if (dataType == "uint16_t")
	{
		return H5T_NATIVE_UINT16;
	}
	else if (dataType == "uint32_t")
	{
		return H5T_NATIVE_UINT32;
	}
	else if (dataType == "uint64_t")
	{
		return H5T_NATIVE_UINT64;
	}
	else
	{
		throw std::runtime_error("Bad data type");
	}
}

// write group attribute
inline void NYXDataLoader::writeGroupAttrib(
	hid_t& in_file_id, std::string const& in_group_name
)
{

	// either default group or not
	auto& group_attrib = groupsAttribs[in_group_name];
	assert(not group_attrib.empty());

	// create a HDF5 group
	hid_t group = H5Gcreate2(
					  in_file_id, ("/" + in_group_name).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
				  );

	auto write = [](hid_t const & attrib, hid_t const & datatype, auto & item)
	{
		// export the given attribute data of arbitrary type
		return (datatype == H5T_NATIVE_INT ? H5Awrite(attrib, datatype, (int*) item.data) :
				(datatype == H5T_NATIVE_DOUBLE ? H5Awrite(attrib, datatype, (double*) item.data) :
				 H5Awrite(attrib, datatype, (char*) item.data)));
	};

	// write each attribute of the group
	for (auto && item : group_attrib)
	{
		hid_t dataspace;
		hid_t attrib;

		// retrieve attribute meta-data
		std::string const& dataname = item.paramName;
		hid_t const datatype = getNativeDataType(item.dataType);
		hsize_t const dims[] = { item.numElements };
		bool const is_array  = (dims[0] > 1);

		// set right dataspace
		if (is_array)
		{
			dataspace = H5Screate(H5S_SIMPLE);
			H5Sset_extent_simple(dataspace, 1, dims, NULL);
			attrib = H5Acreate(group, dataname.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
		}
		else
		{
			assert(datatype == H5T_NATIVE_INT or datatype == H5T_NATIVE_DOUBLE);
			dataspace = H5Screate(H5S_SCALAR);
			attrib = H5Acreate(group, dataname.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT);
		}
		// write it
		herr_t const status = write(attrib, datatype, item);
		log << "\twriting attrib: /" << in_group_name << "/" << item.paramName
			<< ", status: " << (status == -1) << std::endl;

		H5Aclose(attrib);
		H5Sclose(dataspace);
	}

	// clear it to avoid duplicated attribute writing on next compressor
	group_attrib.clear();

	H5Gclose(group);
	log << "Attributes writing finished for: " << in_group_name << std::endl;
}


// no need to pass data_type, as well as group uncompressed data vector
inline void NYXDataLoader::writeGroupData(
	hid_t& in_file_id,
	hid_t& in_filespace,
	hid_t& in_memspace,
	std::string const& in_group_name,
	hsize_t const count[3], hsize_t const offset[3]
)
{

	// should normally be initialized by 'sizePerDim' and 'rankOffset'
	assert(count != nullptr);
	assert(offset != nullptr);

	// step 1: initialization
	// retrieve group fields data
	auto& group_data = (
						   in_group_name == defaultGroup ? toWriteData : groupsData[in_group_name]
					   );

	assert(not group_data.empty());
	int const num_fields = group_data.size();
	std::vector<hid_t> dset_id(num_fields);

	// create a HDF5 group
	hid_t group = H5Gcreate2(
					  in_file_id, ("/" + in_group_name).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
				  );

	// step 2: create metadata
	int field_index = 0;
	for (auto && item : group_data)
	{
		std::string fieldname = "/" + in_group_name + "/" + item.paramName;
		hid_t data_type = getNativeDataType(item.dataType);

		dset_id[field_index] = H5Dcreate(
								   in_file_id, fieldname.c_str(), data_type, in_filespace,
								   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
							   );
		auto const status = (dset_id[field_index] == -1);
		log << "creating fieldname" << fieldname << ", status: " << status << std::endl;
		field_index++;
	}


	// step 3: write fields related data
	hid_t filespace = H5Dget_space(dset_id[0]);
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	field_index = 0;
	for (auto && item : group_data)
	{
		hid_t data_type = getNativeDataType(item.dataType);
		int const shift = count[0] * count[1] * count[2];

		auto write = [&](auto * data)
		{
			herr_t const status = H5Dwrite(
									  dset_id[field_index], data_type, in_memspace, filespace, plist_id, data
								  );

			log << "writing field: " << item.paramName
				<< ", status: " << (status == -1) << std::endl
				<< data[0] << ", " << data[1] << ", "
				<< data[shift - 2] << ", " << data[shift - 1] << std::endl;
		};

		// FIXME: to be refactored, change the way data is stored
		if (item.dataType == "float")
		{
			write((float*) item.data);
		}
		else if (item.dataType == "double")
		{
			write((double*) item.data);
		}
		else if (item.dataType == "int8_t")
		{
			write((int8_t*) item.data);
		}
		else if (item.dataType == "int16_t")
		{
			write((int16_t*) item.data);
		}
		else if (item.dataType == "int32_t")
		{
			write((int32_t*) item.data);
		}
		else if (item.dataType == "int64_t")
		{
			write((int64_t*) item.data);
		}
		else if (item.dataType == "uint8_t")
		{
			write((uint8_t*) item.data);
		}
		else if (item.dataType == "uint16_t")
		{
			write((uint16_t*) item.data);
		}
		else if (item.dataType == "uint32_t")
		{
			write((uint32_t*) item.data);
		}
		else if (item.dataType == "uint64_t")
		{
			write((uint64_t*) item.data);
		}
		else
		{
			log << "Error: Bad data type" << std::endl;
			throw std::runtime_error("Bad data type");
		}
		field_index++;
	}

	// finalization
	herr_t status = H5Gclose(group);
	log << "group close status: " << (status == -1) << std::endl;

	H5Sclose(filespace);
	H5Pclose(plist_id);

	// IMPORTANT
	for (auto && dataset : dset_id)
	{
		H5Dclose(dataset);
	}
	dset_id.clear();

	for (auto && item : group_data)
	{
		item.deAllocateMem();
	}
	group_data.clear();
}

// a new version of 'writeData' that use 'writeGroupData'
inline int NYXDataLoader::writeData(std::string in_filename)
{

	hsize_t fileDims[3];
	hsize_t count[3];
	hsize_t offset[3];

	log << "writing to " << in_filename << std::endl;

	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);

	hid_t file_id = H5Fcreate(in_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

	// write native field data
	// important: we assume here that dimensions and offset per rank
	// are (and should normally be) identical for all datasets (native or derived)
	fileDims[0] = origDims[0];
	fileDims[1] = origDims[1];
	fileDims[2] = origDims[2];

	count[0] = sizePerDim[0];
	count[1] = sizePerDim[1];
	count[2] = sizePerDim[2];

	offset[0] = rankOffset[0];
	offset[1] = rankOffset[1];
	offset[2] = rankOffset[2];

	hid_t filespace = H5Screate_simple(3, fileDims, NULL);
	hid_t memspace  = H5Screate_simple(3, count, NULL);

	log << "count " << count[0] << ", " << count[1] << ", " << count[2] << std::endl;
	log << "offset " << offset[0] << ", " << offset[1] << ", " << offset[2] << std::endl;

	// create attribute 'format'
	{
		char strformat[] = "nyx-lyaf";
		hid_t type = H5Tcopy(H5T_C_S1);
		int32_t len = strlen(strformat);
		H5Tset_size(type, len);
		hid_t ds = H5Screate(H5S_SCALAR);

		hid_t attribute_id = H5Acreate(file_id, "format", type, ds, H5P_DEFAULT, H5P_DEFAULT);
		H5Awrite(attribute_id, type, &strformat);
		H5Aclose(attribute_id);
	}

	// write native fields datasets with the right args
	writeGroupData(file_id, filespace, memspace, defaultGroup, count, offset);

	// write uncompressed derived fields datasets
	for (auto && dataset : groupsData)
	{
		// filespace and memspace should be the same as for native fields
		writeGroupData(file_id, filespace, memspace, dataset.first, count, offset);
	}

	for (auto && group_attrib : groupsAttribs)
	{
		writeGroupAttrib(file_id, group_attrib.first);
	}

	H5Sclose(filespace);
	H5Sclose(memspace);
	H5Pclose(plist_id);
	H5Fclose(file_id);
}

// pass other fields to output HDF5 file
// - process datasets in the first place
// - process single attributes values then
inline bool NYXDataLoader::loadUncompressedFields(nlohmann::json const& jsonInput)
{

	std::string filename = jsonInput["input"]["filename"];
	std::ifstream checkfile(filename);
	if (not checkfile.good())
	{
		log << "\tError: unable to load file: " << filename << std::endl;
		checkfile.close();
		return false;
	}
	checkfile.close();

	hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	/*struct Metadata {
	  std::string name;
	  std::string type;
	  int size;
	};

	// 1. register scalars to be forwarded
	std::unordered_map<std::string, std::vector<std::string>> scalars;
	std::unordered_map<std::string, std::vector<Metadata>> attribs;
	*/

	// parse fields to be stored
	auto const& uncompressed_fields = jsonInput["input"]["uncompressed"];
	for (auto && entry : uncompressed_fields)
	{
		log << "Loading uncompressed group data:" << std::endl;
		assert(entry.find("group") != entry.end());
		std::string const& group = entry["group"];

		if (entry.find("scalars") != entry.end())
		{
			// first handle group datasets
			for (std::string field : entry["scalars"])
			{
				std::string const full_name = "/" + group + "/" + field;
				log << "\tField: " << full_name << std::endl;
				if (loadGroupDataSet(file, group, field) < 0)
				{
					log << "\tError while loading dataset: " << full_name << std::endl;
					return false;
				}
			}
		}
		else if (entry.find("attributes") != entry.end())
		{
			// then handle group attributes
			for (auto && attribute : entry["attributes"])
			{
				/*        attribs[name].emplace_back(
						  attribute["name"], attribute["type"], attribute["size"]
						);*/
				std::string const& name = attribute["name"];
				std::string const& type = attribute["type"];
				int const size = attribute["size"];

				std::string const full_name = "/" + group + "/" + name;
				log << "\tAttribute: " << full_name << std::endl;

				if (loadAttribute(file, group, name, size) < 0)
				{
					log << "\tError while loading attribute: " << full_name << std::endl;
					return false;
				}
			}
		}
	}

	/*
	  // 2. load them
	  log << "Loading group datasets:" << std::endl;
	  for (auto&& entry : scalars) {
		auto const& group = entry.first;
		for (auto&& field : entry.second) {
		  log << "- current field: /" << group << "/" << field << std::endl;
		  loadGroupDataSet(file, group, field);
		}
	  }

	  // 3. load attributes
	  log << "Loading group attributes" << std::endl;
	  for (auto&& entry : attribs) {
		auto const& group = entry.first;
		for (auto&& attrib : entry.second) {
		  log << "- current attribute: /" << group << "/" << attrib << std::endl;
		}
	  }
	*/

	H5Fclose(file);
	//scalars.clear();
	//values.clear();

	// 3. write them
	return true;
}

#endif

/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/

#pragma once

#include <sstream>
#include <iostream>
#include <string>
#include <unordered_map>

#include "dataLoaderInterface.hpp"

#include "strConvert.hpp"
#include "json.hpp"
#include "timer.hpp"
#include "utils.hpp"

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_DataServices.H"  
#include "AMReX_Utility.H"
#include "AMReX_FArrayBox.H"

#include "hdf5.h"

class PltDataLoader: public DataLoaderInterface
{
	int numRanks;
	int myRank;

	std::vector<std::string> uncompressedScalars;
	std::vector<ScalarData> toWriteData;
	bool outputFileCreated;

	amrex::DataServices dataServices;


	int createHDF5(std::string file_path);
	void writeHDF5Metadata(std::string file_path,
		const int nx, const int ny, const int nz, const double domain_size,
		const double omega_b, const double omega_m, const double omega_l,
		const double h, const double z);


	int createHDF5Field(std::string file_path, std::string field_path, std::string units, const int nx, const int ny, const int nz);
	int writeHDF5Field(std::string file_path, std::string field_path, ScalarData d);

	double parseByName(std::string file_name, std::string var_name);

  public:
	PltDataLoader();
	~PltDataLoader();

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int saveCompData(std::string paramName, void * cData);
	int writeData(std::string _filename);

	int saveInputFileParameters() { return 1; }
	int close() { return deAllocateMem(dataType, data); }
	void setParam(std::string paramName, std::string type, std::string value) {};
	bool loadUncompressedFields(nlohmann::json const&) { return false; }



};


inline PltDataLoader::PltDataLoader()
{
	outputFileCreated = false;
}


inline PltDataLoader::~PltDataLoader()
{

}


inline void PltDataLoader::init(std::string _filename, MPI_Comm _comm)
{
	filename = _filename;
	comm = _comm;
	saveData = false;

	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);
	
	dataType = "double";

	amrex::Initialize(comm);	
}



inline int PltDataLoader::loadData(std::string paramName)
{
	log << "loading plt file " << filename << std::endl;
	param = paramName;

	amrex::DataServices::SetBatchMode();
	amrex::Amrvis::FileType fileType(amrex::Amrvis::NEWPLT);
	amrex::DataServices dataServices(filename.c_str(), fileType);
	if (!dataServices.AmrDataOk()) 
	{
		std::cout << "!dataServices.AmrDataOk()" << std::endl;
		amrex::DataServices::Dispatch(amrex::DataServices::ExitRequest, NULL);
	}

	amrex::ParallelDescriptor::IOProcessor();


	//
	// Read metadata.
	//
	amrex::AmrData& amrData = dataServices.AmrDataRef();
	int finestLevel = amrData.FinestLevel();


	int nComp = 1;
	amrex::Vector<int> comps(nComp);
	int i_dm_scalar(amrData.StateNumber(paramName.c_str()));
	comps[0] = i_dm_scalar;
   

	const amrex::Vector<string>& plotVarNames = amrData.PlotVarNames();
	amrex::Vector<string> inVarNames(nComp);
	amrex::Vector<int> destFillComps(nComp);

	amrex::ParallelDescriptor::IOProcessor(); 
  
	for (int i = 0; i < nComp; ++i) 
	{
		inVarNames[i] = plotVarNames[comps[i]];
		amrex::ParallelDescriptor::IOProcessor();
		// if (amrex::ParallelDescriptor::IOProcessor()) 
		// {
		//     std::cout << "    " << amrData.StateNumber(inVarNames[i])
		//               << " (" << inVarNames[i] << ")" << std::endl;
		// }
		destFillComps[i] = i;
	}

	//
	// Make boxes and boxarray.
	//

	// Grab the number of cells in each direction (should be the same, but
	// we're keeping things general).
	amrex::Box pd(amrData.ProbDomain()[0]);
	int64_t grid_nx = pd.bigEnd(0) - pd.smallEnd(0) + 1;
	int64_t grid_ny = pd.bigEnd(1) - pd.smallEnd(1) + 1;
	int64_t grid_nz = pd.bigEnd(2) - pd.smallEnd(2) + 1;
	int64_t num_cells = grid_nx * grid_ny * grid_nz;

	amrex::ParallelDescriptor::IOProcessor();



	 // Check if we can split along x evenly.
	if (grid_nx % numRanks != 0) 
	{
		if (amrex::ParallelDescriptor::IOProcessor()) 
		{
			std::cout << "ERROR: domain decomposition." << std::endl;
			std::cout << "The number of MPI ranks must fit evenly in the number of cells along x." << std::endl;
		}
		MPI_Barrier(comm);
		amrex::Finalize();
	}


	int chunk_size = grid_nx / numRanks;

	// The box for z-pencils
	amrex::Box bx_pencil(amrData.ProbDomain()[0]);

	// List of pencil boxes.
	amrex::BoxList box_list;

	// indexes
	int i, ix_lo, ix_hi;
	for (i=0; i<numRanks; ++i) 
	{
		ix_lo = chunk_size * i;
		ix_hi = ix_lo + chunk_size - 1;

		amrex::Box skewer;
		skewer.setSmall(0, ix_lo);
		skewer.setBig(0, ix_hi);
		skewer.setSmall(1, 0);
		skewer.setBig(1, grid_ny - 1);
		skewer.setSmall(2, 0);
		skewer.setBig(2, grid_nz - 1);

		box_list.push_back(skewer);
	}

	amrex::BoxArray ba(box_list);

	amrex::ParallelDescriptor::IOProcessor(); // Creating multifabs

   

	int num_comps = 1;  // keeping it simple...
	int comp_start = 0;
	int level = 0;
	int ng = 0;         // don't add any ghost zones.

	amrex::DistributionMapping dmap(ba, amrex::ParallelDescriptor::NProcs());

	


	//
	// Start reading.
	//
	amrex::MultiFab mf1(ba, dmap, num_comps, ng);

	amrData.FillVar(mf1, level, paramName.c_str(), comp_start);
	amrData.FlushGrids(amrData.StateNumber(paramName.c_str()));

	if (amrex::ParallelDescriptor::IOProcessor()) 
	{
		for (amrex::MFIter mfi(mf1); mfi.isValid(); ++mfi) 
		{
			// Get fab, lo and hi vectors.
			const amrex::Box& box = mfi.validbox();
			const int *lo_vec = box.loVect();
			const int *hi_vec = box.hiVect();

			// this is in fortran order.
			double *fab_data = mf1[mfi].dataPtr();

			// We must make a copy of the fab data in c-order
			int nx = hi_vec[0] - lo_vec[0] + 1;
			int ny = hi_vec[1] - lo_vec[1] + 1;
			int nz = hi_vec[2] - lo_vec[2] + 1;
			size_t num_cells = nx * ny * nz;


			// Passing information to constructor
			totalNumberOfElements = num_cells;

			dataType = "double";
			elemSize = 8;

			origDims[0] = nx;
			origDims[1] = ny;
			origDims[2] = nz;


			// MPI split!!!
			Partition current = getPartition(myRank, numRanks, origDims[0], origDims[1], origDims[2]);
			rankOffset[0] = current.min_x;
			rankOffset[1] = current.min_y;
			rankOffset[2] = current.min_z;

			sizePerDim[0] = current.max_x - current.min_x;
			sizePerDim[1] = current.max_y - current.min_y;
			sizePerDim[2] = current.max_z - current.min_z;

			numElements = sizePerDim[0]*sizePerDim[1]*sizePerDim[2];


			// Create space for data and store data there
			allocateMem(dataType, numElements, 0, data);


			size_t indexLocal = 0;
			for (size_t z=0; z<sizePerDim[2]; z++)
				for (size_t y=0; y<sizePerDim[1]; y++)
					for (size_t x=0; x<sizePerDim[0]; x++)
					{
						size_t indexGlobal = (rankOffset[2] + z) * (origDims[0] * origDims[1]) +
											 (rankOffset[1] + y) *  origDims[0] +
											 (rankOffset[0] + x);

						// Choose based on datatype
						if (dataType == "float")
							((float *)data)[indexLocal] = fab_data[indexGlobal];
						else if (dataType == "double")
							((double *)data)[indexLocal] = fab_data[indexGlobal];
						else if (dataType == "int")
							((int *)data)[indexLocal] = fab_data[indexGlobal];
						indexLocal++;
					}
		}
	}
	else
	{
		std::cout << myRank << " ~ Something failed" << std::endl;
	}


	log << "original dims: " << origDims[0] << ", " << origDims[1] << ", " << origDims[1] << std::endl;
	log << "dims on this rank: " << sizePerDim[0] << ", " << sizePerDim[1] << ", " << sizePerDim[1] << std::endl;
	log << "rankOffset: " << rankOffset[0] << ", " << rankOffset[1] << ", " << rankOffset[1] << std::endl;
	log << "dataType : " << dataType << std::endl;
	log << "elemSize : " << elemSize << std::endl << std::endl;
	
	log << "loading plt param " << paramName  << " done!"<< std::endl;

	return 1;
}



inline int PltDataLoader::saveCompData(std::string paramName, void * cData)
{
	Timer clock;
	clock.start();
	

	ScalarData temp(paramName, dataType, numElements, sizePerDim, rankOffset);
	temp.allocateMem();
	memcpy(temp.data, cData, numElements*getDataypeSize(dataType));
	toWriteData.push_back(temp);


	clock.stop();
	log << "saving data took: " << clock.getDuration() << " s" << std::endl;
}



inline int PltDataLoader::writeData(std::string _filename)
{
	Timer clock;
	clock.start();

	// Added the extension
	_filename = _filename + ".h5";

	//
	// Create the HDF5 file
	createHDF5(_filename);

	//
	// Add metadata to this file
    amrex::Real comoving_a;
    if (amrex::ParallelDescriptor::IOProcessor()) 
    {
        std::string a_file_path = filename + "/comoving_a";
        std::ifstream a_file;
        a_file.open(a_file_path.c_str(), std::ios::in);
        if (!a_file.good()) 
        {
            amrex::FileOpenFailed(a_file_path);
        }
        a_file >> comoving_a;
    }
	double z = 1.0/comoving_a - 1.0;

	std::string params_file = filename + "/the_parameters";
    double h = parseByName(params_file, "nyx.comoving_h");
    double omega_b = parseByName(params_file, "nyx.comoving_OmB");
    double omega_m = parseByName(params_file, "nyx.comoving_OmM");
    double omega_l = 1.0-omega_m;
    double domain_size = parseByName(params_file, "geometry.prob_hi") * h;

    std::cout << "aaa" << std::endl;

    amrex::DataServices::SetBatchMode();
	amrex::Amrvis::FileType fileType(amrex::Amrvis::NEWPLT);
	amrex::DataServices dataServices(filename, fileType);
	if( ! dataServices.AmrDataOk()) 
	{
    	amrex::DataServices::Dispatch(amrex::DataServices::ExitRequest, NULL);
	}

	amrex::ParallelDescriptor::IOProcessor();

    amrex::AmrData & amrData = dataServices.AmrDataRef();

	amrex::Box pd(amrData.ProbDomain()[0]);
	int64_t grid_nx = pd.bigEnd(0) - pd.smallEnd(0) + 1;
	int64_t grid_ny = pd.bigEnd(1) - pd.smallEnd(1) + 1;
	int64_t grid_nz = pd.bigEnd(2) - pd.smallEnd(2) + 1;
	int64_t num_cells = grid_nx * grid_ny * grid_nz;

	std::cout << "bbb" << std::endl;

	writeHDF5Metadata(_filename.c_str(),
        grid_nx, grid_ny, grid_nz, domain_size,
        omega_b, omega_m, omega_l, h, z);

	

	//
	// Write out the compressed values
	for (int i=0; i<toWriteData.size(); i++)
	{
		std::cout << "ccc " << i << " : " << toWriteData[i].dims[0] << ", " << toWriteData[i].dims[1] << ", " <<toWriteData[i].dims[2] << std::endl;
		createHDF5Field(_filename, toWriteData[i].paramName, "", toWriteData[i].dims[0],toWriteData[i].dims[1],toWriteData[i].dims[2]);
		std::cout << "ddd" << std::endl;

		writeHDF5Field(_filename, toWriteData[i].paramName, toWriteData[i]);

		std::cout << "eee" << std::endl;
	}

	// write out the uncompressed values

	
	clock.stop();
	log << "writing data took: " << clock.getDuration() << " s" << std::endl;
}




inline void PltDataLoader::writeHDF5Metadata(std::string file_path,
	 const int nx, const int ny, const int nz, const double domain_size,
	 const double omega_b, const double omega_m, const double omega_l,
	 const double h, const double z) 
{

	hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

	//
	// Output format root attribute
	//

	const std::string format_text = "nyx-lyaf";

	hid_t str_type     = H5Tcopy(H5T_C_S1);
	hid_t scalar_space = H5Screate(H5S_SCALAR);

	// Fix the str_type length for the format string.
	H5Tset_size(str_type, strlen(format_text.c_str()));

	hid_t attr = H5Acreate(file, "format", str_type, scalar_space, H5P_DEFAULT,
						   H5P_DEFAULT);
	H5Awrite(attr, str_type, format_text.c_str());



	//
	// Domain metadata group
	//

	hid_t group = H5Gcreate(file, "domain", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Dataspace for domain attributes
	hsize_t num_dims[1] = {3};
	hid_t grid_attr_space = H5Screate_simple(1, num_dims, NULL);

	// Grid shape
	int shape[3] = {nx, ny, nz};
	attr = H5Acreate(group, "shape", H5T_STD_I32LE, grid_attr_space,
					 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_NATIVE_INT, shape);

	// Grid size
	double size[3] = {domain_size, domain_size, domain_size};
	attr = H5Acreate(group, "size", H5T_IEEE_F64LE, grid_attr_space,
					 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_NATIVE_DOUBLE, size);


	//
	// Universe metadata group
	//

	group = H5Gcreate(file, "universe", H5P_DEFAULT, H5P_DEFAULT,
					  H5P_DEFAULT);

	attr = H5Acreate(group, "omega_b", H5T_IEEE_F64LE, scalar_space,
					 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_NATIVE_DOUBLE, &omega_b);
	attr = H5Acreate(group, "omega_m", H5T_IEEE_F64LE, scalar_space,
					 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_NATIVE_DOUBLE, &omega_m);
	attr = H5Acreate(group, "omega_l", H5T_IEEE_F64LE, scalar_space,
					 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_NATIVE_DOUBLE, &omega_l);
	attr = H5Acreate(group, "hubble", H5T_IEEE_F64LE, scalar_space,
					 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_NATIVE_DOUBLE, &h);
	attr = H5Acreate(group, "redshift", H5T_IEEE_F64LE, scalar_space,
					 H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_NATIVE_DOUBLE, &z);



	//
	// Field groups
	//

	group = H5Gcreate(file, "native_fields", H5P_DEFAULT, H5P_DEFAULT,
					  H5P_DEFAULT);

	group = H5Gcreate(file, "derived_fields", H5P_DEFAULT, H5P_DEFAULT,
					  H5P_DEFAULT);

	// Close all resources.
	H5Sclose(grid_attr_space);
	H5Gclose(group);
	H5Aclose(attr);
	H5Sclose(scalar_space);
	H5Tclose(str_type);
	H5Fclose(file);
	H5close();
}


inline int PltDataLoader::createHDF5(std::string file_path) 
{
	hid_t file = H5Fcreate(file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file < 0) 
	{
		//std::cout << "Error: could not create file at " << file_path << std::end;
		return -1;
	}

	H5Fclose(file);

	return 1;
}



inline int PltDataLoader::createHDF5Field(std::string file_path, std::string field_path,
	 const std::string units, const int nx, const int ny, const int nz) 
{
	// Open the output.
	hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	// Create a 3D, nx x ny x nz dataspace.
	hsize_t dims[3] = {nx, ny, nz};
	hid_t grid_space = H5Screate_simple(3, dims, NULL);
	// Create the dataset.
	hid_t dataset = H5Dcreate(file, field_path.c_str(), H5T_IEEE_F32LE,
							  grid_space, H5P_DEFAULT, H5P_DEFAULT,
							  H5P_DEFAULT);

	if (dataset < 0) 
	{
		//std::cout << "Error: could not create dataset. H5 returned " << std::end;
		return -1;
	}

	// Don't forget to attach units attr.
	hid_t scalar_space = H5Screate(H5S_SCALAR);
	hid_t str_type = H5Tcopy(H5T_C_S1);
	H5Tset_size(str_type, strlen(units.c_str()));
	hid_t attr = H5Acreate(dataset, "units", str_type, scalar_space,
						   H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, str_type, units.c_str());

	// Close resources.
	H5Aclose(attr);
	H5Tclose(str_type);
	H5Sclose(scalar_space);
	H5Dclose(dataset);
	H5Sclose(grid_space);
	H5Fclose(file);

	return 1;
}


/**
Write the only component in the multifab to the dataset given by field_name.
Uses hdf5-parallel.
*/
inline int PltDataLoader::writeHDF5Field(std::string file_path, std::string field_path, ScalarData d) 
{

	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Info info = MPI_INFO_NULL;
	int mpi_rank;
	MPI_Comm_rank(comm, &mpi_rank);

	// Create the file access prop list.
	hid_t pa_plist = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(pa_plist, comm, info);

	// Open the file, and the group.
	hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDWR, pa_plist);
	// Open the field dataset.
	hid_t dataset = H5Dopen(file, field_path.c_str(), H5P_DEFAULT);

	// Make sure the dataset is there.
	if (dataset < 0) {
		printf("Error on rank %i: Could not find dataset %s.\n", mpi_rank,
			   field_path.c_str());
		exit(1);
	}

	// Grab the dataspace of the field dataset from file.
	hid_t file_dataspace = H5Dget_space(dataset);

	// Create collective io prop list.
	hid_t collective_plist = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(collective_plist, H5FD_MPIO_COLLECTIVE);

	// Iterate over Fabs, select matching hyperslab and write.
	hid_t status;
	// slab lo index and shape.
	hsize_t slab_offsets[3], slab_dims[3];
	hid_t slab_dataspace;

	int write_count = 0;

	//for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
		// Get fab, lo and hi vectors.
		// const Box& box = mfi.validbox();
		// const int *lo_vec = box.loVect();
		// const int *hi_vec = box.hiVect();

		// this is in fortran order.
		//double *fab_data = mf[mfi].dataPtr();

		// We must make a copy of the fab data in c-order
		// int nx = hi_vec[0] - lo_vec[0] + 1;
		// int ny = hi_vec[1] - lo_vec[1] + 1;
		// int nz = hi_vec[2] - lo_vec[2] + 1;

		int nx = d.dims[0];
		int ny = d.dims[1];
		int nz = d.dims[2];
		size_t num_cells = nx * ny * nz;

		// //size_t block_size = num_cells;
		// std::vector<double> h5_data_vec(num_cells);
		// double *h5_data = h5_data_vec.data();
		// //double* h5_data = (double*) malloc(block_size * sizeof(double));
		// if (h5_data == NULL) {
		// 	printf("Error allocating h5 data block.\n");
		// 	exit(1);
		// }

		// //int ix, iy, iz, c_index, f_index;
		// for (int ix = 0; ix < nx; ++ix) {
		// 	for (int iy = 0; iy < ny; ++iy) {
		// 		for (int iz = 0; iz < nz; ++iz) {
		// 			size_t c_index = (ix * ny + iy) * nz + iz;
		// 			size_t f_index = (iz * ny + iy) * nx + ix;
		// 			h5_data[c_index] = fab_data[f_index];
		// 		}
		// 	}
		// }
		double* h5_data = static_cast<double *>(d.data);

		// Data is preped. Now write it.

		// Set slab offset and shape.
		// slab_offsets[0] = lo_vec[0];
		// slab_offsets[1] = lo_vec[1];
		// slab_offsets[2] = lo_vec[2];

		slab_offsets[0] = d.offset[0];
		slab_offsets[1] = d.offset[1];
		slab_offsets[2] = d.offset[2];

		slab_dims[0] = nx;
		slab_dims[1] = ny;
		slab_dims[2] = nz;

		// Create the slab space.
		slab_dataspace = H5Screate_simple(3, slab_dims, NULL);

		// Select the hyperslab matching this fab.
		status = H5Sselect_hyperslab(file_dataspace, H5S_SELECT_SET,
									 slab_offsets, NULL, slab_dims, NULL);
		if (status < 0) {
			printf("Error on rank %i: could not select hyperslab.\n", mpi_rank);
			exit(1);
		}

		// Write this pencil.
		status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, slab_dataspace,
						  file_dataspace, collective_plist, h5_data);
		if (status < 0) {
			printf("Error on rank %i: could not write hyperslab.\n", mpi_rank);
			exit(1);
		}

		H5Sclose(slab_dataspace);
	//	write_count++;
	//}

	// Close HDF5 resources.
	H5Pclose(collective_plist);
	H5Sclose(file_dataspace);
	H5Dclose(dataset);
	H5Fclose(file);
	H5Pclose(pa_plist);
}


inline double PltDataLoader::parseByName(std::string file_name, std::string var_name) 
{
    std::ifstream fin;
    std::string line = "";

    fin.open(file_name.c_str());
    while(std::getline(fin, line)) 
    {
        if (line.find("//") == 0 || line.empty())
            continue;
        else 
        {
            int pos = line.find("=");
            std::string name   = line.substr(0,pos-1);
            std::string svalue = line.substr(pos+1);
            if(name == var_name) 
            {
                std::istringstream instr(svalue);
                double value;
                instr >> value;
                return value;
            }
        }
    }
    return -9999.0;
}
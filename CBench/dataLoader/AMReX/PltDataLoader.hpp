/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/

#pragma once

#include <sstream>
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

class PltDataLoader: public DataLoaderInterface
{
	int numRanks;
	int myRank;

	float origRealDims[3];

	void *tempData;

	amrex::DataServices dataServices;

  public:
	PltDataLoader();
	~PltDataLoader();

	void init(std::string _filename, MPI_Comm _comm);
	int loadData(std::string paramName);
	int saveCompData(std::string paramName, void * cData);
	int writeData(std::string _filename);

	int saveInputFileParameters() { return 1; };
	int close() { return deAllocateMem(dataType, data); }
	void setParam(std::string paramName, std::string type, std::string value) {};
	bool loadUncompressedFields(nlohmann::json const&) { return false; }
};


inline PltDataLoader::PltDataLoader()
{
	amrex::Initialize(MPI_COMM_WORLD);
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


	amrex::DataServices::SetBatchMode();
	amrex::Amrvis::FileType fileType(amrex::Amrvis::NEWPLT);
	dataServices.Init(filename.c_str(), fileType);
	if (!dataServices.AmrDataOk()) 
	{
		std::cout << "!dataServices.AmrDataOk()" << std::endl;
		amrex::DataServices::Dispatch(amrex::DataServices::ExitRequest, NULL);
	}

	if (amrex::ParallelDescriptor::IOProcessor()) 
	{
		std::cout << "done." << std::endl;
	}

	
	dataType = "double";

	

	log << "original dims: " << origDims[0] << ", " << origDims[1] << ", " << origDims[1] << std::endl;
	log << "real dims : " << origRealDims[0] << ", " << origRealDims[1] << ", " << origRealDims[1] << std::endl;
	log << "dataType : " << dataType << std::endl;
}



inline int PltDataLoader::loadData(std::string paramName)
{
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

	// if (amrex::ParallelDescriptor::IOProcessor()) 
	// {
	//     printf("The grid has a shape of (%i, %i, %i), %.3e cells total.\n\n",
	//            (int)grid_nx, (int)grid_ny, (int)grid_nz, (double)num_cells);
	//     printf("Making skewer chunks.\n");
	//     fflush(stdout);
	// }


	 // Check if we can split along x evenly.
	if (grid_nx % myRank != 0) 
	{
		if (amrex::ParallelDescriptor::IOProcessor()) 
		{
			std::cout << "ERROR: domain decomposition." << std::endl;
			std::cout << "The number of MPI ranks must fit evenly in the number of cells along x." << std::endl;
		}
		MPI_Barrier(comm);
		amrex::Finalize();
	}



	int chunk_size = grid_nx / myRank;

	// The box for z-pencils
	amrex::Box bx_pencil(amrData.ProbDomain()[0]);

	// List of pencil boxes.
	amrex::BoxList box_list;

	// indexes
	int i, ix_lo, ix_hi;
	for (i=0; i<myRank; ++i) 
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
	// Start outputting.
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
		std::cout << "Something failed" << std::endl;
	}
}



inline int PltDataLoader::saveCompData(std::string paramName, void * cData)
{
	Timer clock;
	clock.start();
	
	allocateMem(dataType, numElements, 0, tempData);
	memcpy(tempData, cData, numElements*getDataypeSize(dataType));

	clock.stop();
	log << "saving data took: " << clock.getDuration() << " s" << std::endl;
}



inline int PltDataLoader::writeData(std::string _filename)
{
	Timer clock;
	clock.start();

	
	clock.stop();
	log << "writing data took: " << clock.getDuration() << " s" << std::endl;
}



#include <iostream>
#include <vector>

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_DataServices.H"
#include "AMReX_Utility.H"
#include "AMReX_FArrayBox.H"

int main(int argc, char **argv)
{
	int rank, nprocs;
    MPI_Init(&argc,&argv); 
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs); 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); 


    amrex::Initialize(MPI_COMM_WORLD);


    amrex::DataServices::SetBatchMode();
    amrex::Amrvis::FileType fileType(amrex::Amrvis::NEWPLT);
    amrex::DataServices dataServices( "/bigData/plt02669", fileType);
    if (!dataServices.AmrDataOk()) 
    {
    	std::cout << "!dataServices.AmrDataOk()" << std::endl;
        amrex::DataServices::Dispatch(amrex::DataServices::ExitRequest, NULL);
    }


    if (amrex::ParallelDescriptor::IOProcessor()) 
    {
        std::cout << "done." << std::endl;
    }
   


    amrex::AmrData& amrData = dataServices.AmrDataRef();
    int finestLevel = amrData.FinestLevel();

    int Nlev = finestLevel + 1;


    //
    // Read metadata.
    //

    int nComp = 2;
    amrex::Vector<int> comps(nComp);
    //Dark matter fields
    int i_dm_temp(amrData.StateNumber("Temp"));
    comps[0] = i_dm_temp;
    int i_dm_density(amrData.StateNumber("density"));
    comps[1] = i_dm_density;

    amrex::Print() << comps[0];
    amrex::Print() << comps[1];


    const amrex::Vector<string>& plotVarNames = amrData.PlotVarNames();
    amrex::Vector<string> inVarNames(nComp);
    amrex::Vector<int> destFillComps(nComp);
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::cout << std::endl << "Converting the following states: " << std::endl;
    }

    for (int i = 0; i < nComp; ++i) {
        inVarNames[i] = plotVarNames[comps[i]];
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::cout << "    " << amrData.StateNumber(inVarNames[i])
                      << " (" << inVarNames[i] << ")" << std::endl;
        }
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

    if (amrex::ParallelDescriptor::IOProcessor()) {
        printf("The grid has a shape of (%i, %i, %i), %.3e cells total.\n\n",
               (int)grid_nx, (int)grid_ny, (int)grid_nz, (double)num_cells);
        printf("Making skewer chunks.\n");
        fflush(stdout);
    }


     // Check if we can split along x evenly.
    if (grid_nx % nprocs != 0) {
        if (amrex::ParallelDescriptor::IOProcessor()) {
            printf("ERROR: domain decomposition.\n");
            printf("The number of MPI ranks must fit evenly in the number of cells along x.\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        amrex::Finalize();
    }



    int chunk_size = grid_nx / nprocs;

    // The box for z-pencils
    amrex::Box bx_pencil(amrData.ProbDomain()[0]);
    // List of pencil boxes.
    amrex::BoxList box_list;
    // indexes
    int i, ix_lo, ix_hi;
    for (i = 0; i < nprocs; ++i) {
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

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::cout << "Creating multifabs." << std::endl;
    }


    int num_comps = 1;  // keeping it simple...
    int comp_start = 0;
    int level = 0;
    int ng = 0;         // don't add any ghost zones.


    amrex::DistributionMapping dmap(ba, amrex::ParallelDescriptor::NProcs());

    


    //
    // Start outputting.
    //

    // if (amrex::ParallelDescriptor::IOProcessor()) {
    //     std::cout << "Preparing output file." << std::endl;
    // }


    //
    // temperature
    //

    amrex::MultiFab mf1(ba, dmap, num_comps, ng);


    amrData.FillVar(mf1, level, "Temp", comp_start);
    amrData.FlushGrids(amrData.StateNumber("Temp"));

    if (amrex::ParallelDescriptor::IOProcessor()) 
    {
	    std::cout << "valid Temp: " << amrData.StateNumber("Temp") << std::endl;

	    amrex::MFIter mfi(mf1);
	    std::cout << "mfi.isValid(): " << mfi.isValid() << std::endl;


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

	        // this is in fortran order.
	        std::cout << fab_data[0] << std::endl;
	        std::cout << fab_data[1] << std::endl;
	        std::cout << fab_data[2] << std::endl;
	        std::cout << fab_data[3] << std::endl;
	        std::cout << fab_data[4] << std::endl;
	    }

	    // std::string infile; 
	    // std::cout << "Reading " << infile << "...";
	    // pp.get("infile",infile);
	}
	else
	{
		std::cout << "Something failed" << std::endl;
	}



	//
    // density
    //

    amrex::MultiFab mf2(ba, dmap, num_comps, ng);

	amrData.FillVar(mf2, level, "density", comp_start);
    amrData.FlushGrids(amrData.StateNumber("density"));

    if (amrex::ParallelDescriptor::IOProcessor()) 
    {
	    std::cout << "valid Temp: " << amrData.StateNumber("Temp") << std::endl;

	    amrex::MFIter mfi(mf2);
	    std::cout << "mfi.isValid(): " << mfi.isValid() << std::endl;


	    for (amrex::MFIter mfi(mf2); mfi.isValid(); ++mfi) 
	    {
	    	// Get fab, lo and hi vectors.
	        const amrex::Box& box = mfi.validbox();
	        const int *lo_vec = box.loVect();
	        const int *hi_vec = box.hiVect();

	        // this is in fortran order.
	        double *fab_data = mf2[mfi].dataPtr();

	        // We must make a copy of the fab data in c-order
	        int nx = hi_vec[0] - lo_vec[0] + 1;
	        int ny = hi_vec[1] - lo_vec[1] + 1;
	        int nz = hi_vec[2] - lo_vec[2] + 1;
	        size_t num_cells = nx * ny * nz;

	        std::cout << "nx: " << nx << std::endl;
	        std::cout << "ny: " << ny << std::endl;
	        std::cout << "nz: " << nz << std::endl;

	        // this is in fortran order.
	        std::cout << fab_data[0] << std::endl;
	        std::cout << fab_data[1] << std::endl;
	        std::cout << fab_data[2] << std::endl;
	        std::cout << fab_data[3] << std::endl;
	        std::cout << fab_data[4] << std::endl;
	    }

	    // std::string infile; 
	    // std::cout << "Reading " << infile << "...";
	    // pp.get("infile",infile);
	}
	else
	{
		std::cout << "Something failed" << std::endl;
	}

	amrex::ParallelDescriptor::Barrier();
    amrex::Finalize();
    MPI_Finalize(); 
	
	return 0;
}
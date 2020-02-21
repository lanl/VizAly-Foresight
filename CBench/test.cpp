#include <iostream>
#include <vector>

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_DataServices.H"
#include "AMReX_Utility.H"
#include "AMReX_FArrayBox.H"

int main(int argc, char **argv)
{
	// if (argc < 3) {
 //        print_usage(argc, argv);
 //    }

	amrex::Initialize(argc, argv);

    //amrex::ParmParse::Initialize(argc-2, argv+2, "/bigData/plt02669");


    amrex::DataServices::SetBatchMode();
    amrex::Amrvis::FileType fileType(amrex::Amrvis::NEWPLT);
    amrex::DataServices dataServices( "/bigData/plt02669", fileType);

    if (amrex::ParallelDescriptor::IOProcessor()) 
    {
        std::cout << " done." << std::endl;
    }
    else
    {
    	std::cout << "Failed!" << std::endl;
    }


    amrex::AmrData& amrData = dataServices.AmrDataRef();
    int finestLevel = amrData.FinestLevel();
    int Nlev = finestLevel + 1;


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

    int num_comps = 1;  // keeping it simple...
    int comp_start = 0;
    int level = 0;
    int ng = 0;         // don't add any ghost zones.

    amrex::BoxList box_list;


    amrex::BoxArray ba(box_list);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::cout << "Creating multifabs." << std::endl;
    }

    amrex::DistributionMapping dmap(ba, amrex::ParallelDescriptor::NProcs());

    amrex::MultiFab mf1(ba, dmap, num_comps, ng);

    amrData.FillVar(mf1, level, "Temp", comp_start);
    amrData.FlushGrids(amrData.StateNumber("Temp"));

    std::cout << "valid " << amrData.StateNumber("Temp") << std::endl;

    amrex::MFIter mfi(mf1);
    std::cout << "valid " << mfi.isValid() << std::endl;

    for (amrex::MFIter mfi(mf1); mfi.isValid(); ++mfi) 
    {
        // this is in fortran order.
        double *fab_data = mf1[mfi].dataPtr();
        std::cout << fab_data[0] << std::endl;
        std::cout << fab_data[1] << std::endl;
        std::cout << fab_data[2] << std::endl;
        std::cout << fab_data[3] << std::endl;
        std::cout << fab_data[4] << std::endl;
    }

    // std::string infile; 
    // std::cout << "Reading " << infile << "...";
    // pp.get("infile",infile);
	
	return 0;
}
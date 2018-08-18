#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include <mpi.h>

// Helper Functions
#include "json.hpp"
#include "timer.hpp"
#include "log.hpp"
#include "memory.hpp"

// Metrics
#include "relativeError.hpp"
#include "absoluteError.hpp"

// Interfaces
#include "dataLoaderInterface.hpp"
#include "compressorInterface.hpp"

// Readers
#include "thirdparty/genericio/GenericIO.h"
#include "HACCDataLoader.hpp"

// Compressors
#include "blosccompressor.hpp"
#include "BigCrunchcompressor.hpp"


int main(int argc, char *argv[])
{
	// Parse arguments and read json file
	if (argc < 2)
	{
		std::cout << "Input argument needed. Run as: ../inputs/blosc.json" << std::endl;
        std::cout << "Read arguments: " << argc << std::endl;
		return 0;
	}
	

	//
	// init MPI
	int myRank, numRanks, threadSupport;
	MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &threadSupport);
	MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);


	//
	// Read input from json
	nlohmann::json jsonInput;
	std::ifstream jsonFile(argv[1]);
	jsonFile >> jsonInput;


	//
	// Load in the parameters
	std::string inputFileType = jsonInput["input"]["filetype"]; 
	std::string inputFile = jsonInput["input"]["filename"]; 

    std::string outputLogFile = jsonInput["output"]["logfname"];
    outputLogFile += "_" + std::to_string(myRank);

    std::string metricsFile = jsonInput["output"]["metricsfname"];

	std::vector< std::string > params;
	for (int j=0; j<jsonInput["input"]["scalars"].size(); j++)
    	params.push_back( jsonInput["input"]["scalars"][j] );

    std::vector< std::string > compressors;
    for (int j=0; j<jsonInput["compressors"].size(); j++)
        compressors.push_back( jsonInput["compressors"][j] );


    //
    // Create log and metrics files
    std::stringstream debuglog;
    std::stringstream metricslog;

    writeLog(outputLogFile, debuglog.str());

    //
	// Open file
    DataLoaderInterface *ioMgr;

	if (inputFileType == "HACC")
		ioMgr = new HACCDataLoader();
	else
	{
		std::cout << "Unsupported file!!!" << std::endl;
		return 0;
	}

	ioMgr->init(inputFile, MPI_COMM_WORLD);

	std::cout << "Starting compressors...\n";
    // 
    // Cycle through compressors and parameters

    CompressorInterface *compressorMgr;

    // Compressors
    for (int c=0; c<compressors.size(); ++c)
    {

		if (compressors[c] == "blosc")
			compressorMgr = new BLOSCCompressor();
		else if (compressors[c] == "BigCrunch")
			compressorMgr = new BigCrunchCompressor();
		else
        {
            std::cout << "Unsupported compressor: " << compressors[c] << "...Skipping!" << std::endl;
            continue;
        }  

		// init
		compressorMgr->init();

    	// Cycle through params
        for (int i=0; i<params.size(); i++)
        {
            Memory memLoad;
            memLoad.start();

            ioMgr->loadData(params[i]);
            MPI_Barrier(MPI_COMM_WORLD);

            debuglog << ioMgr->getDataInfo();
            debuglog << ioMgr->getLog();
            appendLog(outputLogFile, debuglog);

            MPI_Barrier(MPI_COMM_WORLD);

			//
			// compress
			void * cdata = NULL;
			compressorMgr->compress(ioMgr->data, cdata, ioMgr->type(), ioMgr->size());

			//
			// decompress
			void * decompdata = NULL;
            compressorMgr->decompress(cdata, decompdata, ioMgr->type(), ioMgr->size() );

			writeLogApp(outputLogFile, compressorMgr->getLog());
			compressorMgr->clearLog();

			//
            // Compute metrics
            std::vector<double> rel_err(ioMgr->size());
            std::vector<double> abs_err(ioMgr->size());
                
            for (std::size_t j = 0; j < ioMgr->size(); ++j)
            {
                // Max set tolerence to 1
                double err = relativeError(static_cast<float *>(ioMgr->data)[j], static_cast<float *>(decompdata)[j], 1);
                double err2 = absoluteError(static_cast<float *>(ioMgr->data)[j], static_cast<float *>(decompdata)[j]);
                rel_err.push_back(err);
                abs_err.push_back(err2);
            }

            double max_rel_err = *std::max_element(rel_err.begin(), rel_err.end());
            double max_abs_err = *std::max_element(abs_err.begin(), abs_err.end());

            double total_max_rel_err = 0;
            double total_max_abs_err = 0;
            MPI_Reduce(&max_rel_err, &total_max_rel_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&max_abs_err, &total_max_abs_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            double total_rel_err = total_max_rel_err;
            double total_abs_err = total_max_abs_err;
               
            debuglog << "\n ----- " << params[i] << " error metrics ----- " << std::endl;
            debuglog << " Max Rel Error: " << total_rel_err << std::endl;
            debuglog << " Max Abs Error: " << total_abs_err << std::endl;
            debuglog << "-----------------------------\n";

			std::free(decompdata);

            debuglog << "Memory in use: " << memLoad.getMemoryInUseInMB() << " MB" << std::endl;
            ioMgr->close();

            memLoad.stop();
            debuglog << "Memory leaked: " << memLoad.getMemorySizeInMB() << " MB" << std::endl;

            appendLog(outputLogFile, debuglog);
    		MPI_Barrier(MPI_COMM_WORLD);
    	}

		compressorMgr->close();
    }

    //writeLog(outputLogFile, outputLogFile.str());
    if (myRank == 0)
        std::cout << " Complete! " << std::endl;

	MPI_Finalize();

    

    
	return 0;
}


/*

Run:
mpirun -np 2 CBench ../inputs/jesus_blosc.json

*/

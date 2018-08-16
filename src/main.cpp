#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include <mpi.h>

// Helper Functions
#include "json.hpp"
#include "timer.hpp"
#include "log.hpp"

// Metrics
#include "relativeError.hpp"
#include "absoluteError.hpp"

// Readers
#include "thirdparty/genericio/GenericIO.h"

// Compressors
extern "C"
{
#include <blosc.h>
}

#include "dataLoaderInterface.hpp"
#include "HACCDataLoader.hpp"


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

    std::string inputCompressorType = jsonInput["compressor"];



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

    
    //
	// Cycle through params
    for (int i = 0; i < params.size(); i++)
    {
        ioMgr->loadData(params[i]);
        MPI_Barrier(MPI_COMM_WORLD);

        debuglog << ioMgr->getDataInfo();
        debuglog << ioMgr->getLog();
        appendLog(outputLogFile, debuglog);

        MPI_Barrier(MPI_COMM_WORLD);

        Timer totTime; totTime.start();
        //
        // ---------- Abstract to separate compress/decompress class -----------
        if (inputCompressorType == "blosc")
        {
            //Timer initTime; initTime.start();
            blosc_init();
            //initTime.stop(); metricslog << inputCompressorType << ":InitTime: " << initTime.getDuration() << std::endl;
            // compress
            Timer cTime; cTime.start();
            // Default Input Params: {clevel=9, shuffle=1, sizeof(data), idatasize, idata, cdata, odatasize);
            size_t isize = ioMgr->elemSize*ioMgr->numElements;
            size_t osize = isize + BLOSC_MAX_OVERHEAD;



            void * cdata = std::malloc(isize); //byte array;
            osize = blosc_compress(9, 1, ioMgr->elemSize, isize, ioMgr->data, cdata, osize);
            

            if (osize < 0)
            {
                throw std::runtime_error("Compression error. Error code: " + std::to_string(osize));
            }
            if (osize > 0)
            {
                cdata = std::realloc(cdata, osize);
            }
            debuglog << "\n" << inputCompressorType << " ~ InputBytes: " << isize << ", OutputBytes: " << osize << std::endl;

            cTime.stop(); 
            debuglog << inputCompressorType << " ~ CompressTime: " << cTime.getDuration() << " s " << std::endl;

            // decompress
            Timer dTime; dTime.start();
            osize = ioMgr->elemSize*ioMgr->numElements;
            void * decompdata = std::malloc(osize);
            size_t sz = blosc_decompress(cdata, decompdata, osize);
            if (sz < 0)
                throw std::runtime_error("Decompression error. Error code: " + std::to_string(sz));
            dTime.stop(); debuglog << inputCompressorType << " ~ DeompressTime: " << dTime.getDuration() << " s " << std::endl;

            blosc_destroy();


            // Compute metrics
            std::vector<double> rel_err(ioMgr->numElements);
            std::vector<double> abs_err(ioMgr->numElements);
            
            for (std::size_t j = 0; j < ioMgr->numElements; ++j)
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
        }
        totTime.stop();
        debuglog << "Total Runtime: " << totTime.getDuration() << " s" << std::endl << std::endl << std::endl;

        appendLog(outputLogFile, debuglog);
		MPI_Barrier(MPI_COMM_WORLD);
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

/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
 - Chris Biwer
================================================================================*/

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
	for (int j = 0; j < jsonInput["input"]["scalars"].size(); j++)
		params.push_back( jsonInput["input"]["scalars"][j] );

	std::vector< std::string > compressors;
	for (int j = 0; j < jsonInput["compressors"].size(); j++)
		compressors.push_back( jsonInput["compressors"][j] );


	//
	// Create log and metrics files
	std::stringstream debuglog;
	std::stringstream metricsInfo;

	writeLog(outputLogFile, debuglog.str());

	metricsInfo << "Input file: " << inputFile << std::endl;

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
	// Cycle through compressors and parameters
	CompressorInterface *compressorMgr;

	// Compressors
	for (int c = 0; c < compressors.size(); ++c)
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
		metricsInfo << "\n---------------------------------------" << std::endl;
		metricsInfo << "Compressor: " << compressorMgr->getCompressorName() << std::endl;

		// Cycle through params
		for (int i = 0; i < params.size(); i++)
		{
			Timer compressClock, decompressClock;
			Memory memLoad;

			memLoad.start();

			assert ( ioMgr->loadData(params[i]) == 1);

			MPI_Barrier(MPI_COMM_WORLD);

			debuglog << ioMgr->getDataInfo();
			debuglog << ioMgr->getLog();
			appendLog(outputLogFile, debuglog);

			MPI_Barrier(MPI_COMM_WORLD);

			//
			// compress
			void * cdata = NULL;

			compressClock.start();
			compressorMgr->compress(ioMgr->data, cdata, ioMgr->getTypeSize(), ioMgr->getNumElements());
			compressClock.stop();

			//
			// decompress
			void * decompdata = NULL;

			decompressClock.start();
			compressorMgr->decompress(cdata, decompdata, ioMgr->getTypeSize(), ioMgr->getNumElements() );
			decompressClock.stop();

			writeLogApp(outputLogFile, compressorMgr->getLog());
			compressorMgr->clearLog();


			//
			// Compute metrics
			std::vector<double> rel_err(ioMgr->getNumElements());
			std::vector<double> abs_err(ioMgr->getNumElements());

			for (std::size_t j = 0; j < ioMgr->getNumElements(); ++j)
			{
				// Max set tolerence to 1
				double err = relativeError(static_cast<float *>(ioMgr->data)[j], static_cast<float *>(decompdata)[j], 1);
				double err2 = absoluteError(static_cast<float *>(ioMgr->data)[j], static_cast<float *>(decompdata)[j]);
				rel_err.push_back(err);
				abs_err.push_back(err2);
			}

			double max_rel_err = *std::max_element(rel_err.begin(), rel_err.end());
			double max_abs_err = *std::max_element(abs_err.begin(), abs_err.end());
			double compress_time = compressClock.getDuration();
			double decompress_time = decompressClock.getDuration();
			double compress_throughput = (ioMgr->getNumElements() * ioMgr->getTypeSize()) / compress_time;
			double decompress_throughput = (ioMgr->getNumElements() * ioMgr->getTypeSize()) / decompress_time;

			double total_max_rel_err = 0;
			double total_max_abs_err = 0;
			double max_compress_throughput = 0;
			double max_decompress_throughput = 0;

			MPI_Reduce(&max_rel_err, &total_max_rel_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&max_abs_err, &total_max_abs_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&compress_throughput, &max_compress_throughput, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&decompress_throughput, &max_decompress_throughput, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

			debuglog << "\n ----- " << params[i] << " error metrics ----- " << std::endl;
			debuglog << " Max Rel Error: " << max_rel_err << std::endl;
			debuglog << " Max Abs Error: " << max_abs_err << std::endl;
			debuglog << " Compress time: " << max_abs_err << std::endl;
			debuglog << " Decompress time: " << max_abs_err << std::endl;
			debuglog << "-----------------------------\n";

			metricsInfo << "\nField: " << params[i] << std::endl;
			metricsInfo << "Max Rel Error: " << max_rel_err << std::endl;
			metricsInfo << "Max Abs Error: " << max_abs_err << std::endl;
			metricsInfo << "Compression Throughput: " << compress_throughput << " bytes/s" << std::endl;
			metricsInfo << "DeCompression Throughput: " << decompress_throughput << " bytes/s" << std::endl;

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

	if (myRank == 0)
		writeFile(metricsFile, metricsInfo.str());

	MPI_Finalize();


	return 0;
}


/*

Run:
mpirun -np 2 CBench ../inputs/blosc.json

*/

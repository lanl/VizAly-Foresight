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

// Interfaces
#include "dataLoaderInterface.hpp"
#include "compressorInterface.hpp"
#include "metricInterface.hpp"

// Readers
#include "thirdparty/genericio/GenericIO.h"
#include "HACCDataLoader.hpp"


// Compressors
#include "blosccompressor.hpp"
#ifdef CBENCH_HAS_BIG_CRUNCH
#include "BigCrunchcompressor.hpp"
#endif
#ifdef CBENCH_HAS_SZ
#include "SZcompressor.hpp"
#endif

// Metrics
#include "relativeError.hpp"
#include "absoluteError.hpp"
#include "meansquareError.hpp"

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

	std::vector< std::string > metrics;
	for (int j = 0; j < jsonInput["metrics"].size(); j++)
		metrics.push_back(jsonInput["metrics"][j]);

	//
	// Create log and metrics files
	std::stringstream debuglog;
	std::stringstream metricsInfo;
	std::stringstream csvOutput;
	writeLog(outputLogFile, debuglog.str());

	csvOutput << "Compressor_field" << ", ";
	for (int m = 0; m < metrics.size(); ++m)
		csvOutput << metrics[m] << ", ";
	csvOutput << "Compression Throughput" << ", " << "DeCompression Throughput" << std::endl;
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
	MetricInterface *metricsMgr;

	// Loop compressors
	for (int c = 0; c < compressors.size(); ++c)
	{	
		// Process compressors
		if (compressors[c] == "blosc")
			compressorMgr = new BLOSCCompressor();
	  #ifdef CBENCH_HAS_BIG_CRUNCH
		else if (compressors[c] == "BigCrunch")
			compressorMgr = new BigCrunchCompressor();
	  #endif
	  #ifdef CBENCH_HAS_SZ
		else if (compressors[c] == "SZ")
			compressorMgr = new SZCompressor();
	  #endif
		else
		{
			std::cout << "Unsupported compressor: " << compressors[c] << "...Skipping!" << std::endl;
			continue;
		}


		// Check if the parameters field exist
		if (jsonInput["input"].find("parameters") != jsonInput["input"].end())
		{
			// insert parameter into compressor parameter list
			for (auto it=jsonInput["input"]["parameters"].begin(); it != jsonInput["input"]["parameters"].end(); it++)
				compressorMgr->compressorParameters[it.key()] = strConvert::toStr(it.value());	
		}


		// init
		compressorMgr->init();
		metricsInfo << "\n---------------------------------------" << std::endl;
		metricsInfo << "Compressor: " << compressorMgr->getCompressorName() << std::endl;

		debuglog << "-----------------------------------------" << std::endl;
		debuglog << "Compressor: " << compressorMgr->getCompressorName() << std::endl;

		// Cycle through params
		for (int i=0; i<params.size(); i++)
		{
			Timer compressClock, decompressClock;
			Memory memLoad;

			memLoad.start();

			ioMgr->loadData(params[i]);

			
			debuglog << ioMgr->getDataInfo();
			debuglog << ioMgr->getLog();
			appendLog(outputLogFile, debuglog);
			csvOutput << compressorMgr->getCompressorName() << "_" << ioMgr->getParam() << ", ";
			MPI_Barrier(MPI_COMM_WORLD);


			//
			// compress
			void * cdata = NULL;

			compressClock.start();
			compressorMgr->compress(ioMgr->data, cdata, ioMgr->getType(), ioMgr->getTypeSize(), ioMgr->getNumElements());
			compressClock.stop();

			//
			// decompress
			void * decompdata = NULL;

			decompressClock.start();
			compressorMgr->decompress(cdata, decompdata, ioMgr->getType(), ioMgr->getTypeSize(), ioMgr->getNumElements() );
			decompressClock.stop();

			writeLogApp(outputLogFile, compressorMgr->getLog());
			compressorMgr->clearLog();


			//
			// metrics
			debuglog << "\n----- " << params[i] << " error metrics ----- " << std::endl;
			metricsInfo << "\nField: " << params[i] << std::endl;
			for (int m = 0; m < metrics.size(); ++m)
			{
				if (metrics[m] == "relative_error")
					metricsMgr = new relativeError();
				else if (metrics[m] == "absolute_error")
					metricsMgr = new absoluteError(); 
				else if (metrics[m] == "mse")
					metricsMgr = new meansquareError();
				else
				{
					std::cout << "Unsupported metric: " << metrics[c] << "...Skipping!" << std::endl;
					continue;
				}
				metricsMgr->init(MPI_COMM_WORLD);
				metricsMgr->execute(ioMgr->data, decompdata, ioMgr->getNumElements());
				debuglog << metricsMgr->getLog();
				metricsInfo << metricsMgr->getLog();
				csvOutput << metricsMgr->getGlobalValue() << ", ";
				metricsMgr->close();
			}

			//
			// final timings
			double compress_time = compressClock.getDuration();
			double decompress_time = decompressClock.getDuration();

			debuglog << " Compress time: " << compress_time << std::endl;
			debuglog << " Decompress time: " << decompress_time << std::endl;
			debuglog << "-----------------------------\n";

			double max_compress_throughput = 0;
			double max_decompress_throughput = 0;
			double compress_throughput = ((double)(ioMgr->getNumElements() * ioMgr->getTypeSize()) / (1024.0*1024.0) )/ compress_time;
			double decompress_throughput = ((double)(ioMgr->getNumElements() * ioMgr->getTypeSize()) / (1024.0*1024.0) )/ decompress_time;
			MPI_Reduce(&compress_throughput, &max_compress_throughput, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&decompress_throughput, &max_decompress_throughput, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

			metricsInfo << "Compression Throughput: " << max_compress_throughput << " Mbytes/s" << std::endl;
			metricsInfo << "DeCompression Throughput: " << max_decompress_throughput << " Mbytes/s" << std::endl;

			csvOutput << max_compress_throughput << ", " << max_decompress_throughput << "\n";

			//
			// deallocate
			std::free(decompdata);

			debuglog << "Memory in use: " << memLoad.getMemoryInUseInMB() << " MB" << std::endl;
			ioMgr->close();

			memLoad.stop();
			debuglog << "Memory leaked: " << memLoad.getMemorySizeInMB() << " MB \n\n" << std::endl;

			appendLog(outputLogFile, debuglog);
			MPI_Barrier(MPI_COMM_WORLD);
		}

		compressorMgr->close();
	}

	if (myRank == 0)
	{
		writeFile(metricsFile, metricsInfo.str());
		writeFile(metricsFile + ".csv", csvOutput.str());
	}

	MPI_Finalize();


	return 0;
}


/*

Run:
mpirun -np 2 CBench ../inputs/all.json

*/

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

// Factories
#include "compressorFactory.hpp"
#include "metricFactory.hpp"

// Readers
#include "BinaryDataLoader.hpp"
#include "HACCDataLoader.hpp"
#ifdef CBENCH_HAS_NYX
	#include "NYXDataLoader.hpp"
#endif


int main(int argc, char *argv[])
{
	// Parse arguments and read json file
	if (argc < 2)
	{
		std::cout << "Input argument needed. Run as: ../inputs/all.json" << std::endl;
		std::cout << "Read arguments: " << argc << std::endl;
		return 0;
	}

	//
	// init MPI
	int myRank, numRanks, threadSupport;
	MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &threadSupport);
	MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	// For humans
	if (myRank == 0)
		std::cout << "Starting ... look at the log for progress update ... " << std::endl;


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

	csvOutput << "Max Compression Throughput, Max DeCompression Throughput, ";
	csvOutput << "Min Compression Throughput, Min DeCompression Throughput, Compression Ratio" << std::endl;
	metricsInfo << "Input file: " << inputFile << std::endl;



	//
	// Open file
	DataLoaderInterface *ioMgr;

	if (inputFileType == "Binary")
		ioMgr = new BinaryDataLoader();
	else if (inputFileType == "HACC")
		ioMgr = new HACCDataLoader();
  #ifdef CBENCH_HAS_NYX
	else if (inputFileType == "NYX")
		ioMgr = new NYXDataLoader();
  #endif
	else
	{
		if (myRank == 0)
			std::cout << "Unsupported format " << inputFileType << "!!!" << std::endl;
		return 0;
	}

	// Check if the datainfo field exist for a dataset
	if (jsonInput["input"].find("datainfo") != jsonInput["input"].end())
	{
		// insert datainfo into loader parameter list
		for (auto it = jsonInput["input"]["datainfo"].begin(); it != jsonInput["input"]["datainfo"].end(); it++)
			ioMgr->loaderParams[it.key()] = strConvert::toStr(it.value());
	}

	ioMgr->init(inputFile, MPI_COMM_WORLD);

	

	//
	// Cycle through compressors and parameters
	CompressorInterface *compressorMgr;
	MetricInterface *metricsMgr;

	// Loop compressors
	for (int c = 0; c < compressors.size(); ++c)
	{	
		// initialize compressor
		compressorMgr = CompressorFactory::createCompressor(compressors[c]);
		if (compressorMgr == NULL)
		{
			if (myRank == 0)
				std::cout << "Unsupported compressor: " << compressors[c] << " ... Skipping!" << std::endl;
			continue;
		}

		// Check if the parameters field exist
		if (jsonInput["input"].find("parameters") != jsonInput["input"].end())
		{
			// insert parameter into compressor parameter list
			for (auto it=jsonInput["input"]["parameters"].begin(); it != jsonInput["input"]["parameters"].end(); it++)
				compressorMgr->compressorParameters[it.key()] = strConvert::toStr(it.value());	
		}


		// log
		compressorMgr->init();
		metricsInfo << "\n---------------------------------------" << std::endl;
		metricsInfo << "Compressor: " << compressorMgr->getCompressorName() << std::endl;

		debuglog << "===============================================" << std::endl;
		debuglog << "Compressor: " << compressorMgr->getCompressorName() << std::endl;


		// Cycle through params
		for (int i=0; i<params.size(); i++)
		{
			Timer compressClock, decompressClock;
			Memory memLoad;

			memLoad.start();

			ioMgr->loadData(params[i]);

			
			// log stuff
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


			// Get compression ratio
			unsigned long totalCompressedSize;
			unsigned long compressedSize = (unsigned long)compressorMgr->getCompressedSize();
			MPI_Allreduce(&compressedSize, &totalCompressedSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
			
			unsigned long totalUnCompressedSize;
			unsigned long unCompressedSize = ioMgr->getTypeSize() * ioMgr->getNumElements();
			MPI_Allreduce(&unCompressedSize, &totalUnCompressedSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);


			debuglog << "\n\ncompressedSize: " << compressedSize << ", totalCompressedSize: " << totalCompressedSize << std::endl;
			debuglog << "unCompressedSize: " << unCompressedSize << ", totalUnCompressedSize: " << totalUnCompressedSize << std::endl;
			debuglog << "Compression ratio: " << totalUnCompressedSize/(float)totalCompressedSize << std::endl;

			writeLogApp(outputLogFile, compressorMgr->getLog());
			compressorMgr->clearLog();


			//
			// metrics
			debuglog << "\n----- " << params[i] << " error metrics ----- " << std::endl;
			metricsInfo << "\nField: " << params[i] << std::endl;
			for (int m = 0; m < metrics.size(); ++m)
			{
				metricsMgr = MetricsFactory::createMetric(metrics[m]);
				if (metricsMgr == NULL)
				{
					if (myRank == 0)
						std::cout << "Unsupported metric: " << metrics[m] << " ... Skipping!" << std::endl;
					continue;
				}

				metricsMgr->init(MPI_COMM_WORLD);
				metricsMgr->execute(ioMgr->data, decompdata, ioMgr->getNumElements());

				debuglog << metricsMgr->getLog();
				metricsInfo << metricsMgr->getLog();
				csvOutput << metricsMgr->getGlobalValue() << ", ";

				metricsMgr->close();
			}
			debuglog << "-----------------------------\n";
			debuglog << "\nMemory in use: " << memLoad.getMemoryInUseInMB() << " MB" << std::endl;


			//
			// Metrics Computation
			double compress_time = compressClock.getDuration();
			double decompress_time = decompressClock.getDuration();

			double compress_throughput = ((double)(ioMgr->getNumElements() * ioMgr->getTypeSize()) / (1024.0*1024.0) )/ compress_time;
			double decompress_throughput = ((double)(ioMgr->getNumElements() * ioMgr->getTypeSize()) / (1024.0*1024.0) )/ decompress_time;


			double max_compress_throughput = 0;
			double max_decompress_throughput = 0;
			double min_compress_throughput = 0;
			double min_decompress_throughput = 0;

			MPI_Reduce(&compress_throughput, &max_compress_throughput, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
			MPI_Reduce(&decompress_throughput, &max_decompress_throughput, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

			MPI_Reduce(&compress_throughput, &min_compress_throughput, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
			MPI_Reduce(&decompress_throughput, &min_decompress_throughput, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		

			//
			// deallocate
			std::free(decompdata);
			ioMgr->close();
			memLoad.stop();

		
			//
			// log stuff
			debuglog << "\nCompress time: " << compress_time << std::endl;
			debuglog << "Decompress time: " << decompress_time << std::endl;
			debuglog << "\nMemory leaked: " << memLoad.getMemorySizeInMB() << " MB" << std::endl;
			debuglog << "........................................." << std::endl << std::endl;

			appendLog(outputLogFile, debuglog);
		
			if (myRank == 0)
			{
				metricsInfo << "Max Compression Throughput: " << max_compress_throughput << " Mbytes/s" << std::endl;
				metricsInfo << "Max DeCompression Throughput: " << max_decompress_throughput << " Mbytes/s" << std::endl;
				metricsInfo << "Min Compression Throughput: " << min_compress_throughput << " Mbytes/s" << std::endl;
				metricsInfo << "Min DeCompression Throughput: " << min_decompress_throughput << " Mbytes/s" << std::endl;
				metricsInfo << "Compression ratio: " << totalUnCompressedSize/(float)totalCompressedSize << std::endl;
				csvOutput << max_compress_throughput << ", " << max_decompress_throughput << ", " 
						  << min_compress_throughput << ", " << min_decompress_throughput << ", " 
						  << totalUnCompressedSize/(float)totalCompressedSize << "\n";

				writeFile(metricsFile, metricsInfo.str());
				writeFile(metricsFile + ".csv", csvOutput.str());
			}
		
			MPI_Barrier(MPI_COMM_WORLD);
		}
		compressorMgr->close();
	}

	// for humans
	if (myRank == 0)
		std::cout << "That's all folks!" << std::endl;

	MPI_Finalize();

	return 0;
}


/*
Run:
mpirun -np 2 CBench ../inputs/HACC_all.json
*/

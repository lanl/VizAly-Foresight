/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
 - Chris Biwer
==========================================================================*/
#include <iostream>
#include <sstream>
#include <vector>
#include <time.h>
#include <stdlib.h>

#include <mpi.h>

// Helper Functions
#include "json.hpp"
#include "timer.hpp"
#include "log.hpp"
//#include "debugLog.hpp"
#include "memory.hpp"
#include "utils.hpp"

// Interfaces
#include "dataLoaderInterface.hpp"
#include "compressorInterface.hpp"
#include "metricInterface.hpp"

// Factories
#include "compressorFactory.hpp"
#include "metricFactory.hpp"
#include "dataLoaderFactory.hpp"

//std::stringstream Log::Logging::debugLog;


int main(int argc, char *argv[])
{
	//
	// init MPI
	int myRank, numRanks, threadSupport;
	MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &threadSupport);
	MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);


	//
	// Validate input params
	if (!validateInput(argc, argv, myRank, numRanks))
	{
		MPI_Finalize();
		return 0;
	}

	//
	// Load input


	//
	// Pass JSON file to json parser and
	nlohmann::json jsonInput;
	std::ifstream jsonFile(argv[1]);
	jsonFile >> jsonInput;


	//
	// Load in the global input parameters
	std::string inputFilename = jsonInput["input"]["filename"];


	// timesteps
	int minTimestep = 0;
	int maxTimestep = 1;
	if (jsonInput["input"].contains("timesteps"))
	{
		minTimestep = jsonInput["input"]["timesteps"][0];  	// inclusive
		maxTimestep = jsonInput["input"]["timesteps"][1];	// exclusive
	}
	int numTimesteps = maxTimestep - minTimestep;


	// write out decompressed files + output name
	bool writeData = false;
	std::string outputFilename = "";
	if (jsonInput["data-reduction"]["cbench-output"].contains("output-decompressed"))
	{
		writeData = jsonInput["data-reduction"]["cbench-output"]["output-decompressed"];
		if (writeData)
			outputFilename = extractFileName(inputFilename);

		// Initial a random number in case output name is not provided
		srand(time(NULL));
	}


	// location of decompressed files
	std::string outputPath = ".";
	if ( jsonInput["data-reduction"]["cbench-output"].contains("output-decompressed-location"))
	{
		outputPath = jsonInput["data-reduction"]["cbench-output"]["output-decompressed-location"];
		if (myRank == 0)
			createFolder(outputPath);
	}


	// Log file name
	createFolder("logs");
	std::string outputLogFilename = "logs/" + jsonInput["data-reduction"]["cbench-output"]["log-file"].get<std::string>() + "_" + std::to_string(myRank);


	// Store compressors, scalars, and metrics to use
	std::vector<std::string> compressors;
	for (int i = 0; i < jsonInput["data-reduction"]["cbench-compressors"].size(); i++)
		compressors.push_back(jsonInput["data-reduction"]["cbench-compressors"][i]["name"]);

	std::vector<std::string> scalars;
	for (int i = 0; i < jsonInput["input"]["scalars"].size(); i++)
		scalars.push_back(jsonInput["input"]["scalars"][i]);

	std::vector<std::string> metrics;
	for (int i = 0; i < jsonInput["data-reduction"]["cbench-metrics"].size(); i++)
		metrics.push_back(jsonInput["data-reduction"]["cbench-metrics"][i]["name"]);


	//
	// For humans; all seems valid, let's start ...
	if (myRank == 0)
		std::cout << "Starting ... \nLook at the log for progress update ... \n" << std::endl;



	//
	// Create log and metrics files
	Timer overallClock;
	std::stringstream debuglog, metricsInfo, csvOutputHeader;

	overallClock.start();

	
	for (int ts=minTimestep; ts<maxTimestep; ts++)
	{
		//
		// Create header for metrics output
		std::stringstream csvOutput;
		csvOutput << "Compressor_field" << "__" << "params" << ", " << "name, ";
		for (int m = 0; m < metrics.size(); ++m)
			csvOutput << metrics[m] << ", ";
		csvOutput << "Compression Throughput(MB/s), DeCompression Throughput(MB/s), Compression Ratio" << std::endl;

		debuglog << "\n********************************** timestep: " << ts << " of " << numTimesteps << std::endl;


		//
		// Open data file
		std::string fileToLoad;
		DataLoaderInterface *ioMgr;
		{
			//
			// Get filename
			fileToLoad = inputFilename;
			if (numTimesteps > 1)
			{
				std::string tempStr = inputFilename;
				fileToLoad = tempStr.replace( tempStr.find("%"), tempStr.find("%")+1, strConvert::toStr(ts) );
			}
				
			metricsInfo << "Input file: " << fileToLoad << std::endl;
			if (myRank == 0)
				std::cout << "\nReading " << fileToLoad << std::endl;
			

			//
			// Create Data Loader
			std::string inputFileType = jsonInput["input"]["filetype"];
			ioMgr = DataLoaderFactory::createLoader(inputFileType);

			if (ioMgr == NULL)
			{
				if (myRank == 0)
					std::cout << "Unsupported loader: " << inputFileType << " ... exiting!" << std::endl;
				
				debuglog << "Unsupported loader: " << inputFileType << " ... exiting!"<< std::endl;
				writeFile(outputLogFilename, debuglog.str());

				MPI_Finalize();
				return 0;
			}
			else if (inputFileType == "NYX")
			{
				if (jsonInput["input"].contains("group"))
					ioMgr->setParam("group", "string", jsonInput["input"]["group"]);
			}


			//
			// Check if the datainfo field exist for a dataset
			if (jsonInput["input"].contains("datainfo"))
			{
				// insert datainfo into loader parameter list
				for (auto it = jsonInput["input"]["datainfo"].begin(); it != jsonInput["input"]["datainfo"].end(); it++)
					ioMgr->loaderParams[it.key()] = strConvert::toStr(it.value());
			}

			ioMgr->init(fileToLoad, MPI_COMM_WORLD);
			ioMgr->setTimestep(ts);
			ioMgr->setSave(writeData);


			// Save parameters of input file to facilitate rewrite
			if (writeData)
				ioMgr->saveInputFileParameters();
		}


		//
		// Cycle through compressors
		CompressorInterface *compressorMgr;
		for (int c = 0; c < compressors.size(); ++c)
		{
			// initialize compressor
			compressorMgr = CompressorFactory::createCompressor(compressors[c]);
			if (compressorMgr == NULL)
			{
				if (myRank == 0)
					std::cout << "Unsupported compressor: " << compressors[c] << " ... skipping!" << std::endl;

				debuglog << "Unsupported compressor: " << compressors[c] << " ... skipping!" << std::endl;

				continue;
			}


			// initialize compressor
			compressorMgr->init();


			// Apply parameter if same for all scalars, else delay for later
			bool sameCompressorParams = true;
			if (jsonInput["data-reduction"]["cbench-compressors"][c].find("compressor-params") != jsonInput["data-reduction"]["cbench-compressors"][c].end())
				sameCompressorParams = false;
			else
			{
				for (auto it = jsonInput["data-reduction"]["cbench-compressors"][c].begin(); it != jsonInput["data-reduction"]["cbench-compressors"][c].end(); ++it)
					if ((it.key() != "name") && (it.key() != "output-prefix"))
						compressorMgr->compressorParameters[it.key()] = strConvert::toStr(it.value());
			}


			// log
			metricsInfo << "\n---------------------------------------" << std::endl;
			metricsInfo << "Compressor: " << compressorMgr->getCompressorName() << std::endl;

			debuglog << "===============================================" << std::endl;
			debuglog << "Compressor: " << compressorMgr->getCompressorName() << std::endl;


			//
			// Cycle through scalars
			for (int i = 0; i < scalars.size(); i++)
			{
				Timer compressClock, decompressClock;
				Memory memLoad(true);

				// Check if parameter is valid before proceding
				if (!ioMgr->loadData(scalars[i]))
				{
					memLoad.stop();
					std::cout << "ioMgr->loadData(" << scalars[i] << ") failed!" << std::endl; 
					debuglog << "ioMgr->loadData(" << scalars[i] << ") failed!" << std::endl; 
					continue;
				}


				// Read in compressor parameter for this field
				if (!sameCompressorParams)
				{
					compressorMgr->compressorParameters.clear();  // reset compression param for each field
					int numdifferentParams = jsonInput["data-reduction"]["cbench-compressors"][c]["compressor-params"].size();

					for (int cp = 0; cp < numdifferentParams; cp++)
					{
						for (auto it = jsonInput["data-reduction"]["cbench-compressors"][c]["compressor-params"][cp]["scalar"].begin();
								it != jsonInput["data-reduction"]["cbench-compressors"][c]["compressor-params"][cp]["scalar"].end(); it++)
						{
							if (*it != scalars[i])
								continue;

							for (auto itt = jsonInput["data-reduction"]["cbench-compressors"][c]["compressor-params"][cp].begin();
									itt != jsonInput["data-reduction"]["cbench-compressors"][c]["compressor-params"][cp].end(); ++itt)
								if (itt.key() != "scalar")
									compressorMgr->compressorParameters[itt.key()] = strConvert::toStr(itt.value());
						}
					}
				}


				// log stuff
				debuglog << ioMgr->getDataInfo();
				debuglog << ioMgr->getLog();
				ioMgr->clearLog();
				//writeLog(outputLogFilename, debuglog);
				writeFile(outputLogFilename, debuglog.str());

				metricsInfo << compressorMgr->getParamsInfo() << std::endl;
				csvOutput << compressorMgr->getCompressorName() << "_" << scalars[i] << "__" << compressorMgr->getParamsInfo()
						<< ", " << jsonInput["data-reduction"]["cbench-compressors"][c]["output-prefix"].get<std::string>() << ", ";

				MPI_Barrier(MPI_COMM_WORLD);


				//
				// compress
				void *cdata = NULL;

				compressClock.start();
				compressorMgr->compress(ioMgr->data, cdata, ioMgr->getType(), ioMgr->getTypeSize(), ioMgr->getSizePerDim());
				compressClock.stop();


				//
				// decompress
				void *decompdata = NULL;

				decompressClock.start();
				compressorMgr->decompress(cdata, decompdata, ioMgr->getType(), ioMgr->getTypeSize(), ioMgr->getSizePerDim());
				decompressClock.stop();


				// Get compression ratio
				unsigned long totalCompressedSize;
				unsigned long compressedSize = (unsigned long) compressorMgr->getCompressedSize();
				MPI_Allreduce(&compressedSize, &totalCompressedSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

				unsigned long totalUnCompressedSize;
				unsigned long unCompressedSize = ioMgr->getTypeSize() * ioMgr->getNumElements();
				MPI_Allreduce(&unCompressedSize, &totalUnCompressedSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);


				debuglog << "\n\ncompressedSize: " << compressedSize << ", totalCompressedSize: " << totalCompressedSize << std::endl;
				debuglog << "unCompressedSize: " << unCompressedSize << ", totalUnCompressedSize: " << totalUnCompressedSize << std::endl;
				debuglog << "Compression ratio: " << totalUnCompressedSize / (float) totalCompressedSize << std::endl;

				//appendLog(outputLogFilename, compressorMgr->getLog());
				debuglog << compressorMgr->getLog();
				compressorMgr->clearLog();


				//
				// Cycle through metrics
				debuglog << "\n----- " << scalars[i] << " error metrics ----- " << std::endl;
				metricsInfo << "\nField: " << scalars[i] << std::endl;

				MetricInterface *metricsMgr;
				for (int m=0; m<metrics.size(); ++m)
				{
					metricsMgr = MetricsFactory::createMetric(metrics[m]);
					if (metricsMgr == NULL)
					{
						if (myRank == 0)
							std::cout << "Unsupported metric: " << metrics[m] << " ... Skipping!" << std::endl;

						debuglog << "Unsupported metric: " << metrics[m] << " ... Skipping!"<< std::endl;
						continue;
					}


					// Read in additional params for metrics
					for (auto it = jsonInput["data-reduction"]["cbench-metrics"][m].begin(); it != jsonInput["data-reduction"]["cbench-metrics"][m].end(); it++)
					{
						if (it.key() != "name")
							for (auto it2 = jsonInput["data-reduction"]["cbench-metrics"][m][it.key()].begin();
									it2 != jsonInput["data-reduction"]["cbench-metrics"][m][it.key()].end(); it2++)
								if (*it2 != scalars[i])
									continue;
								else
								{
									metricsMgr->parameters[it.key()] = strConvert::toStr(scalars[i]);
									break;
								}
					}


					// Launch
					metricsMgr->init(MPI_COMM_WORLD);
					metricsMgr->execute(ioMgr->data, decompdata, ioMgr->getNumElements());
					debuglog << metricsMgr->getLog();
					metricsInfo << metricsMgr->getLog();
					metricsMgr->clearLog();
					csvOutput << metricsMgr->getGlobalValue() << ", ";

					// darw histogram if needed
					if (myRank == 0)
						if (metricsMgr->additionalOutput != "")
						{
							createFolder("logs");
							std::string outputHistogramName = "logs/";
							outputHistogramName += extractFileName(fileToLoad) + "_" + compressors[c] + "_" + scalars[i];
							outputHistogramName += "_" + metrics[m] + "_" + compressorMgr->getParamsInfo() + "_hist.py";
							writeFile(outputHistogramName, metricsMgr->additionalOutput);
						}
					metricsMgr->close();
				} // metrics
				debuglog << "-----------------------------\n";
				debuglog << "\nMemory in use: " << memLoad.getMemoryInUseInMB() << " MB" << std::endl;



				//
				// Metrics Computation
				double compress_time = compressClock.getDuration();
				double decompress_time = decompressClock.getDuration();

				double compress_throughput   = ((double) (ioMgr->getNumElements() * ioMgr->getTypeSize()) / (1024.0 * 1024.0)) / compress_time;     // MB/s
				double decompress_throughput = ((double) (ioMgr->getNumElements() * ioMgr->getTypeSize()) / (1024.0 * 1024.0)) / decompress_time;	// MB/s


				double max_compress_throughput = 0;
				double max_decompress_throughput = 0;
				double min_compress_throughput = 0;
				double min_decompress_throughput = 0;
				double max_compress_time = 0;

				MPI_Reduce(&compress_throughput, &max_compress_throughput, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
				MPI_Reduce(&decompress_throughput, &max_decompress_throughput, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
				MPI_Reduce(&compress_time, &max_compress_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

				MPI_Reduce(&compress_throughput, &min_compress_throughput, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
				MPI_Reduce(&decompress_throughput, &min_decompress_throughput, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);


				//
				// Store the uncompressed data pointer, does not actually write yet
				if (writeData)
				{
					debuglog << "writing: " << scalars[i] << std::endl;

					ioMgr->saveCompData(scalars[i], decompdata);
					
					debuglog << ioMgr->getLog();
					ioMgr->clearLog();
				}


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
				debuglog << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl << std::endl;

				//writeLog(outputLogFilename, debuglog);
				writeFile(outputLogFilename, debuglog.str());


				//
				// log to metrics file and output metrics csv
				if (myRank == 0)
				{
					metricsInfo << "Max Compression Throughput: " << max_compress_throughput << " Mbytes/s" << std::endl;
					metricsInfo << "Max DeCompression Throughput: " << max_decompress_throughput << " Mbytes/s" << std::endl;
					metricsInfo << "Max Compress Time: " << max_compress_time << " s" << std::endl;

					metricsInfo << "Min Compression Throughput: " << min_compress_throughput << " Mbytes/s" << std::endl;
					metricsInfo << "Min DeCompression Throughput: " << min_decompress_throughput << " Mbytes/s" << std::endl;

					metricsInfo << "Compression ratio: " << totalUnCompressedSize / (float) totalCompressedSize << std::endl;
					csvOutput << min_compress_throughput << ", " << min_decompress_throughput << ", " 
							<< totalUnCompressedSize / (float) totalCompressedSize << "\n";

					{
						std::string metricsFile = jsonInput["data-reduction"]["cbench-output"]["metrics-file"];
						writeFile(metricsFile, metricsInfo.str());

						if (numTimesteps > 1)
							writeFile(metricsFile + std::to_string(ts) + ".csv", csvOutput.str());
						else
							writeFile(metricsFile + ".csv", csvOutput.str());
					}
				}

				MPI_Barrier(MPI_COMM_WORLD);
			}  // scalars


			//
			// write data to disk if requested in the json file
			if (writeData)
			{
				Timer clockWrite;
				clockWrite.start();

				ioMgr->loadUncompressedFields(jsonInput);

				//   #if CBENCH_HAS_NYX
				// 	debuglog << "\nLoading uncompressed fields" << std::endl;
				// 	ioMgr->loadUncompressedFields(jsonInput);
				// 	//debuglog << ioMgr->getLog();
				// 	MPI_Barrier(MPI_COMM_WORLD);
				//   #endif

				// debuglog << "Write data .... \n";

				// // Pass through original data to preserve original file data structure
				// for (int i = 0; i < ioMgr->inOutData.size(); i++)
				// {
				// 	if (!ioMgr->inOutData[i].doWrite)
				// 	{
				// 		debuglog << "writing uncoompressed" << std::endl;
				// 		ioMgr->loadData(ioMgr->inOutData[i].name);
				// 		ioMgr->saveCompData(ioMgr->inOutData[i].name, ioMgr->data);
				// 		ioMgr->close();
				// 	}
				// }

				// std::string decompressedOutputName;
				// if (jsonInput["data-reduction"]["cbench-compressors"][c].find("output-prefix") != jsonInput["data-reduction"]["cbench-compressors"][c].end())
				// 	decompressedOutputName = jsonInput["data-reduction"]["cbench-compressors"][c]["output-prefix"];
				// else
				// 	decompressedOutputName = "__" + compressorMgr->getCompressorName() + "_" + std::to_string(rand());

				// if (outputPath != ".")
				// 	decompressedOutputName = outputPath + "/" + decompressedOutputName + "__" + outputFilename;
				// else
				// 	decompressedOutputName = decompressedOutputName + "__" + outputFilename;

				// // Write out uncompressed (lossy) data
				// ioMgr->writeData(decompressedOutputName);


				//
				// Get the name and path of the new file
				std::string decompressedOutputName, fileToOutput;
				if (jsonInput["compressors"][c].contains("output-prefix"))
					decompressedOutputName = jsonInput["compressors"][c]["output-prefix"];
				else
					decompressedOutputName = "__" + compressorMgr->getCompressorName() + "_" + std::to_string(rand());

				// deal with timesteps
				if (outputFilename.find("%") != std::string::npos)
				{
					std::string tempStr = outputFilename;
					fileToOutput = tempStr.replace( outputFilename.find("%"), outputFilename.find("%")+1, strConvert::toStr(ts) );
				}
				else
					fileToOutput = outputFilename;

				// path for folders
				if (outputPath != ".")
					decompressedOutputName = outputPath + "/" + decompressedOutputName + "__" + fileToOutput;
				else
					decompressedOutputName = decompressedOutputName + "__" + fileToOutput;


				//
				// Write out uncompressed (lossy) data
				ioMgr->writeData(decompressedOutputName);


				clockWrite.stop();

				if (myRank == 0)
					std::cout << "wrote out " << decompressedOutputName << "." << std::endl;

				debuglog << ioMgr->getLog();	ioMgr->clearLog();

				debuglog << "Write output took: " << clockWrite.getDuration() << " s " << std::endl;
				//writeLog(outputLogFilename, debuglog);
				writeFile(outputLogFilename, debuglog.str());
			}  // write Data

			compressorMgr->close();
		} // compressors
	} // timesteps

	overallClock.stop();
	debuglog << "\nTotal run time: " << overallClock.getDuration() << " s " << std::endl;
	//writeLog(outputLogFilename, debuglog);
	writeFile(outputLogFilename, debuglog.str());


	// For humans
	if (myRank == 0)
		std::cout << "\nThat's all folks!" << std::endl;

	MPI_Finalize();

	return 0;
}


/*
Run:
mpirun -np 4 CBench ../inputs/hacc/HACC_all.json
mpirun -np 4 CBench ../inputs/nyx/nyx_all.json
*/




/*=========================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset

==========================================================================*/
#include <iostream>
#include <sstream>
#include <vector>
#include <memory>
#include <time.h>
#include <stdlib.h>

#include <mpi.h>

// Helper Functions
#include "json.hpp"
#include "timer.hpp"
#include "log.hpp"
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


// Global log
std::stringstream debugLog;



//
// Parses the JSON file for the main inputs
inline int parseIput(nlohmann::json jsonInput, 
			int myRank,
			std::string &inputFilename,
			int &numTimesteps,
			std::vector<std::string> &filenameTs,
			int &minTimestep,
			int &maxTimestep,
			bool &writeData,
			std::string &outputLogFilename,
			std::vector<std::string> &compressors,
			std::vector<std::string> &scalars,
			std::vector<std::string> &metrics)
{
	//
	// Load in the global input parameters


	// multi timesteps: TODO: Not fully supported yet!
	inputFilename = "";
	minTimestep = 0;
	maxTimestep = 1;
	if (jsonInput["input"].contains("timesteps"))	// files in order
	{
		minTimestep = jsonInput["input"]["timesteps"][0];  	// inclusive
		maxTimestep = jsonInput["input"]["timesteps"][1];	// exclusive

		inputFilename = jsonInput["input"]["filename"];
	}
	else if (jsonInput["input"].contains("filename-timesteps"))	// arbitrary file names
		{
			maxTimestep = jsonInput["input"]["filename-timesteps"].size();
			for (int i=0; i<maxTimestep; i++)
				filenameTs.push_back( jsonInput["input"]["filename-timesteps"][i] );
		}
		else
			inputFilename = jsonInput["input"]["filename"];	// single timestep

	numTimesteps = maxTimestep - minTimestep;
	if (numTimesteps < 1)
	{
		std::cout << "No timesteps found!!!" << std::endl;
		return 0;
	}

	if (inputFilename == "")
	{
		std::cout << "No input file found!!!" << std::endl;
		return 0;
	}



	// write out decompressed files + output name
	writeData = false;
	if (jsonInput["data-reduction"]["cbench-output"].contains("output-decompressed"))
		writeData = jsonInput["data-reduction"]["cbench-output"]["output-decompressed"];


	// Log file name
	outputLogFilename = "logs/" + jsonInput["data-reduction"]["cbench-output"]["log-file"].get<std::string>() + "_" + std::to_string(myRank);


	// Store compressors, scalars, and metrics to use
	for (int i = 0; i < jsonInput["data-reduction"]["cbench-compressors"].size(); i++)
		compressors.push_back(jsonInput["data-reduction"]["cbench-compressors"][i]["name"]);

	for (int i = 0; i < jsonInput["input"]["scalars"].size(); i++)
		scalars.push_back(jsonInput["input"]["scalars"][i]);

	for (int i = 0; i < jsonInput["data-reduction"]["cbench-metrics"].size(); i++)
		metrics.push_back(jsonInput["data-reduction"]["cbench-metrics"][i]["name"]);


	return 1;
}



inline int initDataLoader(DataLoaderInterface* &ioMgr,
							nlohmann::json jsonInput,
							int myRank,
							int ts,
							std::string &fileToLoad,
							std::string inputFilename,
							int numTimesteps,
							std::vector<std::string> filenameTs,
							std::stringstream &metricsInfo)
{
	//
	// Get filename
	fileToLoad = inputFilename;
	if (numTimesteps > 1)
		if (filenameTs.size() > 0)
			fileToLoad = filenameTs[ts];
		else
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
		
		debugLog << "Unsupported loader: " << inputFileType << " ... exiting!"<< std::endl;
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
	
	return 1;
}




inline void writeDecompressedData(nlohmann::json jsonInput,
									int myRank,
									DataLoaderInterface* ioMgr,
									CompressorInterface* compressorMgr,
									std::string fileToLoad,
									int index)
{	
	Timer clock;
	clock.start("write");

	if (myRank == 0)
		std::cout << "Writing data ... " << std::endl;


	ioMgr->loadUncompressedFields(jsonInput);


	// Pass through original data to preserve original file data structure
	for (int i = 0; i < ioMgr->inOutData.size(); i++)
	{
		if (!ioMgr->inOutData[i].doWrite)
		{
			debugLog << "writing uncoompressed" << std::endl;
			ioMgr->loadData(ioMgr->inOutData[i].name);
			ioMgr->saveCompData(ioMgr->inOutData[i].name, ioMgr->data);
			ioMgr->close();
		}
	}



	//
	// Get the name and path of the new file
	std::string decompressedOutputName, fileToOutput;
	if (jsonInput["data-reduction"]["cbench-compressors"][index].contains("output-prefix"))
		decompressedOutputName = jsonInput["data-reduction"]["cbench-compressors"][index]["output-prefix"];
	else
		decompressedOutputName = "__" + compressorMgr->getCompressorName() + "_" + std::to_string(rand());

	fileToOutput = extractFileName(fileToLoad);				

	// path for folders
	std::string outputPath = ".";
	if ( jsonInput["data-reduction"]["cbench-output"].contains("output-decompressed-location"))
	{
		createFolder(outputPath);
		outputPath = jsonInput["data-reduction"]["cbench-output"]["output-decompressed-location"];
		decompressedOutputName = outputPath + "/" + decompressedOutputName + "__" + fileToOutput;
		//std::cout << outputPath << std::endl;
	}
	else
		decompressedOutputName = decompressedOutputName + "__" + fileToOutput;
	
	//
	// Write out uncompressed (lossy) data
	ioMgr->writeData(decompressedOutputName);


	clock.stop("write");

	if (myRank == 0)
		std::cout << "wrote out " << decompressedOutputName << "." << std::endl;

	debugLog << "Write output took: " << clock.getDuration("write") << " s " << std::endl;
	
}


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


	// Initial a random number generator for filenames
	srand(time(NULL));


	// Status of returned functions: O=bad, 1=good!
	int status = 1;


	//
	// Load inputs
	//

	//
	// Pass JSON file to json parser
	nlohmann::json jsonInput;
	std::ifstream jsonFile(argv[1]);
	jsonFile >> jsonInput;

	//
	// Load in the global input parameters
	std::string inputFilename = "";
	int numTimesteps = 1;
	std::vector<std::string> filenameTs;
	int minTimestep = 0;
	int maxTimestep = 1;
	bool writeData = false;
	std::string outputLogFilename = "log";
	std::vector<std::string> compressors;
	std::vector<std::string> scalars;
	std::vector<std::string> metrics;

	status = parseIput(jsonInput,
			myRank, 
			inputFilename,
			numTimesteps,
			filenameTs,
			minTimestep,
			maxTimestep,
			writeData,
			outputLogFilename,
			compressors,
			scalars,
			metrics);

	writeLog(outputLogFilename, debugLog.str());
	
	assertStatus(status);


	// Creating folders
	if (myRank == 0)
	{
		createFolder("logs");
		std::cout << "Starting ... \nLook at the log for progress update ... \n" << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);		// Syncing on all inputs

	//
	// All seems valid, let's start ...
	//

	//
	// Create log and metrics files
	Timer clock;
	std::stringstream metricsInfo, csvOutputHeader;
	
	clock.start("overall");

	
	for (int ts=minTimestep; ts<maxTimestep; ts++)
	{
		//
		// Create header for metrics output
		std::stringstream csvOutput;
		csvOutput << "Compressor_field" << "__" << "params" << ", " << "name, ";
		for (int m = 0; m < metrics.size(); ++m)
			csvOutput << metrics[m] << ", ";

		csvOutput << "Compression Throughput(MB/s), DeCompression Throughput(MB/s), Compression Ratio" << std::endl;
		debugLog << "\n****************************** timestep: " << ts << " of " << numTimesteps << std::endl;


		//
		// initialize data loader
		std::string fileToLoad;
		DataLoaderInterface *ioMgr;

		status = initDataLoader(ioMgr,
						  		jsonInput,
								myRank,
								ts,
								fileToLoad,
								inputFilename,
								numTimesteps,
								filenameTs,
								metricsInfo);

		writeLog(outputLogFilename, debugLog.str());
		assertStatus(status);

		if (myRank == 0)
			std::cout << "\nReading " << fileToLoad  << " done!" << std::endl;
		MPI_Barrier(MPI_COMM_WORLD);	// All reading synchronized!!!
			



		//
		// Cycle through compressors
		CompressorInterface *compressorMgr;
		for (int c = 0; c < compressors.size(); ++c)
		{
			//
			// initialize compressor
			compressorMgr = CompressorFactory::createCompressor(compressors[c]);
			if (compressorMgr == NULL)
			{
				if (myRank == 0)
					std::cout << "Unsupported compressor: " << compressors[c] << " ... skipping!" << std::endl;

				debugLog << "Unsupported compressor: " << compressors[c] << " ... skipping!" << std::endl;
				continue;
			}
			compressorMgr->init();


			//
			// Apply parameter if same for all scalars, else delay for later
			bool sameCompressorParams = true;
			if (jsonInput["data-reduction"]["cbench-compressors"][c].contains("compressor-params"))
			{
				sameCompressorParams = false;
				debugLog << "sameCompressorParams = false" << std::endl;
			}
			else
			{
				for (auto it  = jsonInput["data-reduction"]["cbench-compressors"][c].begin(); 
						  it != jsonInput["data-reduction"]["cbench-compressors"][c].end(); ++it)
					if ((it.key() != "name") && (it.key() != "output-prefix"))
					{
						compressorMgr->compressorParameters[it.key()] = strConvert::toStr(it.value());
					}
				
				debugLog << "sameCompressorParams = true" << std::endl;
			}


			if (myRank == 0)
				std::cout << "\nCompressor " << compressors[c] << " initialized!" << std::endl;


			// log
			metricsInfo << "\n---------------------------------------" << std::endl;
			metricsInfo << "Compressor: " << compressorMgr->getCompressorName() << std::endl;

			debugLog << "===============================================" << std::endl;
			debugLog << "Compressor: " << compressorMgr->getCompressorName() << std::endl;

			writeLog(outputLogFilename, debugLog.str());

			//
			// Cycle through scalars
			for (int i = 0; i < scalars.size(); i++)
			{
				MPI_Barrier(MPI_COMM_WORLD);
				if (myRank == 0)
					std::cout << "\nLoading " << scalars[i] << " ... " << std::endl;

				Memory memLoad(true);

				// Check if parameter is valid before proceding
				if (!ioMgr->loadData(scalars[i]))
				{
					memLoad.stop();
					std::cout << "ioMgr->loadData(" << scalars[i] << ") failed!" << std::endl; 
					debugLog  << "ioMgr->loadData(" << scalars[i] << ") failed!" << std::endl; 
					writeLog(outputLogFilename, debugLog.str());
					continue;
				}

				// Read in compressor parameter for this field if not the same
				bool scalarFound = false;
				if (!sameCompressorParams)
				{
					compressorMgr->compressorParameters.clear();  // reset compression param for each field
					int numDifferentParams = jsonInput["data-reduction"]["cbench-compressors"][c]["compressor-params"].size();

					for (int cp = 0; cp < numDifferentParams; cp++)
					{
						for (auto it  = jsonInput["data-reduction"]["cbench-compressors"][c]["compressor-params"][cp]["scalar"].begin();
								  it != jsonInput["data-reduction"]["cbench-compressors"][c]["compressor-params"][cp]["scalar"].end(); it++)
						{
							if (*it != scalars[i])
								continue;

							scalarFound = true;
							for (auto itt = jsonInput["data-reduction"]["cbench-compressors"][c]["compressor-params"][cp].begin();
									 itt != jsonInput["data-reduction"]["cbench-compressors"][c]["compressor-params"][cp].end(); ++itt)
								if (itt.key() != "scalar")
								{
									compressorMgr->compressorParameters[itt.key()] = strConvert::toStr(itt.value());
								}
						}
					}
				}

				// log stuff
				debugLog << ioMgr->getDataInfo();
				writeLog(outputLogFilename, debugLog.str());

				MPI_Barrier(MPI_COMM_WORLD); // Sync after reading



				// jinzhenw 06/15/2022
				//std::cout << "scalar name: " << scalars[i] << std::endl;
				float  fmin_o = std::numeric_limits<float>::max();
				float  fmax_o = -1 * std::numeric_limits<float>::max();
				float fmean_o = 0;
				int nrow =  ioMgr->getSizePerDim()[0];
				int ncol =  ioMgr->getSizePerDim()[1];
				int ndep =  ioMgr->getSizePerDim()[2];

				float* data = (float*)ioMgr->data;
				minMeanMax(data, fmin_o, fmean_o, fmax_o, nrow, ncol, ndep);
				//std::cout << "fmin: " << fmin << ", fmean: " << fmean << ", fmax: " << fmax << std::endl;
				

				// Do Compression
				void *decompdata = NULL;
				if (!sameCompressorParams && scalarFound == false)
				{
					// NO compression will be applied here!!!
					if (myRank == 0)
						std::cout << "No Compression!" << std::endl;


					metricsInfo << "None" << std::endl;
					csvOutput   << "None" << "_" << scalars[i] << "__" << "None"
							    << ", " << jsonInput["data-reduction"]["cbench-compressors"][c]["output-prefix"].get<std::string>() << ", ";

					clock.start("compress");
					clock.start("decompress");

					clock.stop("compress");
					clock.stop("decompress");



					size_t numel = ioMgr->getSizePerDim()[0];
					for (int i=1; i<5; i++)
						if (ioMgr->getSizePerDim()[i] != 0)
							numel *= ioMgr->getSizePerDim()[i];

					decompdata = malloc(numel*ioMgr->getTypeSize());
					memcpy(decompdata, ioMgr->data, numel*ioMgr->getTypeSize() );
					compressorMgr->setCompressedSize(numel*ioMgr->getTypeSize());
				}
				else
				{
					if (myRank == 0)
					{
						std::cout << "Compressing " << scalars[i] << "... " << std::endl;
						// jinzhenw 06/20/22
						//std::cout << compressorMgr->getParamsInfo() << std::endl;
					}

					metricsInfo << compressorMgr->getParamsInfo() << std::endl;
					csvOutput << compressorMgr->getCompressorName() << "_" << scalars[i] << "__" << compressorMgr->getParamsInfo()
							<< ", " << jsonInput["data-reduction"]["cbench-compressors"][c]["output-prefix"].get<std::string>() << ", ";

					
					//
					// compress
					void *cdata = NULL;

					clock.start("compress");
					compressorMgr->compress(ioMgr->data, cdata, ioMgr->getType(), ioMgr->getTypeSize(), ioMgr->getSizePerDim());
					clock.stop("compress");


					//
					// decompress

					clock.start("decompress");
					compressorMgr->decompress(cdata, decompdata, ioMgr->getType(), ioMgr->getTypeSize(), ioMgr->getSizePerDim());
					clock.stop("decompress");

				}
				
				MPI_Barrier(MPI_COMM_WORLD);	// sync after compression

				//jinzhenw 06/15/2022
				float  fmin_c = std::numeric_limits<float>::max();
				float  fmax_c = -1 * std::numeric_limits<float>::max();
				float fmean_c = 0;
				data = (float*)decompdata;
				minMeanMax(data, fmin_c, fmean_c, fmax_c, nrow, ncol, ndep);
				//std::cout << "fmin: " << fmin << ", fmean: " << fmean << ", fmax: " << fmax << std::endl;
				//
				float Fmin_o, Fmax_o, Fmean_o;
				float Fmin_c, Fmax_c, Fmean_c;
				MPI_Allreduce(&fmin_o,  &Fmin_o, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
				MPI_Allreduce(&fmax_o,  &Fmax_o, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
				MPI_Allreduce(&fmean_o, &Fmean_o, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
				MPI_Allreduce(&fmin_c,  &Fmin_c, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
				MPI_Allreduce(&fmax_c,  &Fmax_c, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
				MPI_Allreduce(&fmean_c, &Fmean_c, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

				writeLog(outputLogFilename, debugLog.str());
				MPI_Barrier(MPI_COMM_WORLD);	// sync after compression

				//
				// Computing Metrics
				//

				if (myRank == 0)
				{
					std::cout << "Computing metrics ... " << std::endl;
					std::cout << "fmin_o: " << Fmin_o << ", fmean_o: " << Fmean_o/numRanks << ", fmax_o: " << Fmax_o << std::endl;
					std::cout << "fmin_c: " << Fmin_c << ", fmean_c: " << Fmean_c/numRanks << ", fmax_c: " << Fmax_c << std::endl;
				}

				// Get compression ratio
				unsigned long totalCompressedSize;
				unsigned long compressedSize = (unsigned long) compressorMgr->getCompressedSize();
				MPI_Allreduce(&compressedSize, &totalCompressedSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

				unsigned long totalUnCompressedSize;
				unsigned long unCompressedSize = ioMgr->getTypeSize() * ioMgr->getNumElements();
				MPI_Allreduce(&unCompressedSize, &totalUnCompressedSize, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);


				debugLog << "\n\ncompressedSize: " << compressedSize << ", totalCompressedSize: " << totalCompressedSize << std::endl;
				debugLog << "unCompressedSize: "   << unCompressedSize << ", totalUnCompressedSize: " << totalUnCompressedSize << std::endl;
				debugLog << "Compression ratio: "  << totalUnCompressedSize / (float) totalCompressedSize << std::endl;


				//
				// Cycle through metrics
				debugLog << "\n----- " << scalars[i] << " error metrics ----- " << std::endl;
				metricsInfo << "\nField: " << scalars[i] << std::endl;

				MetricInterface *metricsMgr;
				for (int m=0; m<metrics.size(); ++m)
				{
					metricsMgr = MetricsFactory::createMetric(metrics[m]);
					if (metricsMgr == NULL)
					{
						if (myRank == 0)
							std::cout << "Unsupported metric: " << metrics[m] << " ... Skipping!" << std::endl;

						debugLog << "Unsupported metric: " << metrics[m] << " ... Skipping!"<< std::endl;
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
					metricsMgr->execute(ioMgr->data, decompdata, ioMgr->getNumElements(), ioMgr->getType());
					
					csvOutput << metricsMgr->getGlobalValue() << ", ";

					// draw histogram if needed
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
				debugLog << "-----------------------------\n";
				debugLog << "\nMemory in use: " << memLoad.getMemoryInUseInMB() << " MB" << std::endl;


				//
				// Metrics Computation
				double compress_time = clock.getDuration("compress");
				double decompress_time = clock.getDuration("decompress");

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
					debugLog << "writing: " << scalars[i] << std::endl;
					ioMgr->saveCompData(scalars[i], decompdata);
				}


				//
				// deallocate
				std::free(decompdata);


				ioMgr->close();
				memLoad.stop();


				//
				// log stuff
				debugLog << "\nCompress time: " << compress_time << std::endl;
				debugLog << "Decompress time: " << decompress_time << std::endl;
				debugLog << "\nMemory leaked: " << memLoad.getMemorySizeInMB() << " MB" << std::endl;
				debugLog << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl << std::endl;

				writeLog(outputLogFilename, debugLog.str());


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
				if (myRank == 0)
					std::cout << scalars[i] << " processing done!" << std::endl;
				MPI_Barrier(MPI_COMM_WORLD);

			}  // scalars


			//
			// write data to disk if requested in the json file
			if (writeData)
				writeDecompressedData(jsonInput,
					 	  				myRank,
					 	  				ioMgr,
						  				compressorMgr,
						  				fileToLoad,
						  				c);

			writeLog(outputLogFilename, debugLog.str());


			compressorMgr->close();

	
			if (myRank == 0)
				std::cout << compressors[c] << " done!" << std::endl;
		} // compressors
	} // timesteps

	//overallClock.stop();
	clock.stop("overall");
	debugLog << "\nTotal run time: " << clock.getDuration("overall") << " s " << std::endl;
	writeLog(outputLogFilename, debugLog.str());


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




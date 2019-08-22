/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/
#pragma once

#include <stdio.h> 
#include <stdbool.h> 
#include <string>
#include <fstream>
#include <iostream>
#include <sys/stat.h>

#ifdef _WIN32
	#include <direct.h>
	#define mkdir(a, b) _mkdir(a)
#endif

#include "json.hpp"



struct stat finfo;


inline int createFolder(std::string folderName)
{
  #ifdef _WIN32
	folderName = ".\\" + folderName;
  #endif

	int res = stat(folderName.c_str(), &finfo);

	if (finfo.st_mode & S_IFDIR)
		return 1; // Directory already exists 
	else if (res != 0)
	{
		const int dir_err = mkdir(folderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (-1 == dir_err)
		{
			std::cout << "Could not create directory " << folderName << std::endl;
			return 0;
		}
	}
	else 
	{
		std::cout << folderName << " is not a directory!" << std::endl;
		return 0;
	}
	return 1;
}



inline bool fileExisits(char *filename) 
{
	std::ifstream ifs(filename);
	return ifs.good();
}



inline bool isPowerOfTwo(int n) 
{ 
	if (n == 0) 
		return 0; 

	while (n != 1) 
	{ 
		if (n%2 != 0) 
			return 0; 
		n = n/2; 
	} 
	return 1; 
} 



inline std::string extractFileName(std::string inputString)
{
	std::size_t pos = inputString.find_last_of("/\\");
	return inputString.substr(pos+1); 
}



inline int allocateMem(std::string dataType, size_t numElements, int offset, void *& data)
{
	if (dataType == "float")
		data = new float[numElements + offset];
	else if (dataType == "double")
		data = new double[numElements + offset];
	else if (dataType == "int8_t")
		data = new int8_t[numElements + offset];
	else if (dataType == "int16_t")
		data = new int16_t[numElements + offset];
	else if (dataType == "int32_t")
		data = new int32_t[numElements + offset];
	else if (dataType == "int64_t")
		data = new int64_t[numElements + offset];
	else if (dataType == "uint8_t")
		data = new uint8_t[numElements + offset];
	else if (dataType == "uint16_t")
		data = new uint16_t[numElements + offset];
	else if (dataType == "uint32_t")
		data = new uint32_t[numElements + offset];
	else if (dataType == "uint64_t")
		data = new uint64_t[numElements + offset];
	else
		return 0;

	return 1;
}



inline int deAllocateMem(std::string dataType, void *& data)
{
	if (data == NULL) // already deallocated!
		return 1;

	if (dataType == "float")
		delete[](float*) data;
	else if (dataType == "double")
		delete[](double*) data;
	else if (dataType == "int8_t")
		delete[](int8_t*) data;
	else if (dataType == "int16_t")
		delete[](int16_t*) data;
	else if (dataType == "int32_t")
		delete[](int32_t*) data;
	else if (dataType == "int64_t")
		delete[](int64_t*) data;
	else if (dataType == "uint8_t")
		delete[](uint8_t*) data;
	else if (dataType == "uint16_t")
		delete[](uint16_t*) data;
	else if (dataType == "uint32_t")
		delete[](uint32_t*) data;
	else if (dataType == "uint64_t")
		delete[](uint64_t*) data;
	else
		return 0;

	data = NULL;

	return 1;
}




// Function to validate whatever could be wrong with the input
inline int validateInput(int argc, char *argv[], int myRank, int numRanks)
{
	// Check if we have the right number of arguments
	if (argc < 2)
	{
		if (myRank == 0)
		{
			std::cerr << "Input argument needed. Run as: <path to input JSON file>" << std::endl;
			std::cerr << "Read arguments: " << argc << std::endl;
		}

		return 0;
	}

	// Check if input file provided exists
	if (!fileExisits(argv[1]))
	{
		if (myRank == 0)
			std::cerr << "Could not find input JSON file " << argv[1] << "." << std::endl;

		return 0;
	}



	// Validate JSON file
	nlohmann::json jsonInput;
	std::ifstream jsonFile(argv[1]);

	try
	{
		jsonFile >> jsonInput;
	}
	catch (nlohmann::json::parse_error &e)
	{
		if (myRank == 0)
			std::cerr << "Input JSON file " << argv[1] << " is invalid!\n"
					  << e.what() << "\n"
					  << "Validate your JSON file using e.g. https://jsonformatter.curiousconcept.com/ " << std::endl;
		return 0;
	}



	// Check if powers of 2 number of ranks
	if (jsonInput["cbench"]["output"].find("output-decompressed-location") != jsonInput["cbench"]["output"].end())
		if (!isPowerOfTwo(numRanks))
		{
			if (myRank == 0)
				std::cerr << "Please run with powers of 2 number of ranks. e.g. 4, 8, 16, 32, ... when writing out HACC files."
						  << std::endl;

			return 0;
		}

	return 1;
}

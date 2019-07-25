/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/
#pragma once

#include <dataLoaderInterface.hpp>
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


inline bool fileExists(char *filename) 
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

inline int writeCompressedStream(void* input, size_t cbytes, std::string outfile, std::string compressorName, std::string dataType, size_t dataTypeSize, size_t* n, std::stringstream &log)
{
	// Header:
	// compressor-name,ioMgr->getType(), ioMgr->getTypeSize(), ioMgr->getSizePerDim() (5 dims)
	// Data:
	// raw-bytes

	return 1;
}

inline int writeCompressedStream(void* input, size_t cbytes, std::string outfile, std::string compressorName, DataLoaderInterface *ioMgr, std::stringstream&log)
{
	// Header:
	// compressor-name, cbytes, ioMgr->getTypeSize(), ioMgr->getSizePerDim() (5 dims)
	// [c][c][c],[size_t],[size_t],([size_t][size_t][size_t][size_t][size_t]),[bytes]
	// Data:
	// raw-bytes
	size_t dataTypeSize = ioMgr->getTypeSize();
	size_t* n = ioMgr->getSizePerDim();

	std::ofstream out;

	log << "Creating file " << outfile << std::endl;
	out.open(outfile, std::ifstream::binary);

	char* tmp_c = new char[1];
	//tmp_c[0] = ' '; tmp_c[1] = ' '; tmp_c[2] = ' ';
	for (size_t i = 0; i < 3; i++)
	{
		if (i >= compressorName.size())
			tmp_c[0] = ' ';
		else
			tmp_c[0] = compressorName.c_str()[i];
		out.write(reinterpret_cast<char*>(tmp_c), sizeof(char));
	}
	// Write header information
	out.write(reinterpret_cast<char*>(&cbytes), sizeof(size_t));
	out.write(reinterpret_cast<char*>(&dataTypeSize), sizeof(size_t));
	out.write(reinterpret_cast<char*>(&n[0]), sizeof(size_t));
	out.write(reinterpret_cast<char*>(&n[1]), sizeof(size_t));
	out.write(reinterpret_cast<char*>(&n[2]), sizeof(size_t));
	out.write(reinterpret_cast<char*>(&n[3]), sizeof(size_t));
	out.write(reinterpret_cast<char*>(&n[4]), sizeof(size_t));
	
	log << "Wrote: " << 3 + 8 + 8 + (8 * 5) << " byte header";
	// Write Data
	out.write(reinterpret_cast<char*>(input), cbytes);
	out.close();

	log << " and " << cbytes << " bytes of data" << std::endl;

	return 1;
}

inline int readCompressedStream(std::string infile, void* &output, size_t cbytes, DataLoaderInterface* ioMgr, std::stringstream& log)
{
	log << "Read: " << 3 + 8 + 8 + (8 * 5) << " byte header";
	// Write Data

	log << " and " << cbytes << " bytes of data" << std::endl;

	return 1;
}
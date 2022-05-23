/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/

#pragma once

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>


extern std::stringstream debugLog;	// Global log for all files


inline void writeFile(std::string filename, std::string log)
{
	std::ofstream outputFile( filename.c_str(), std::ios::out);
	outputFile << log;
	outputFile.close();
}

inline void writeLog(std::string filename, std::string log)
{
	std::ofstream outputFile( (filename + ".log").c_str(), std::ios::out);
	outputFile << log;
	outputFile.close();
}

/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/

#ifndef _LOG_H_
#define _LOG_H_

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>


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

inline void writeLogApp(std::string filename, std::string log)
{
	std::ofstream outputFile( (filename + ".log").c_str(), std::ios::out | std::ios::app);
	outputFile << log;
	outputFile.close();
}



inline void writeLog(std::string filename, std::stringstream log)
{
	std::ofstream outputFile( (filename + ".log").c_str(), std::ios::out);
	outputFile << log.str();
	outputFile.close();
}

inline void appendLog(std::string filename, std::stringstream & log)
{
	std::ofstream outputFile( (filename + ".log").c_str(), std::ios::out | std::ios::app);
	outputFile << log.str();
	outputFile.close();

	log.str("");	// clears the log
}

#endif

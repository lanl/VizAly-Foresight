#pragma once

#include <sstream>
#include <iostream>

namespace log
{
	static std::stringstream debugLog;

	inline void writeLog(std::string filename, std::stringstream log)
	{
		std::ofstream outputFile( (filename + ".log").c_str(), std::ios::out);
		outputFile << log.str();
		outputFile.close();
	}

	inline void clearLog(std::string filename, std::stringstream log)
	{
		log.str("");
		log.clear();
	}
}

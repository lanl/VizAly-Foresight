#ifndef _LOG_H_
#define _LOG_H_

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>

class Log
{
	std::string outputFilename;
	std::string logMsg;

  public:
	Log();
	Log(std::string _outputFilename): outputFilename(_outputFilename){ logMsg = ""; }
	~Log();

	void setOutputFilename(std::string _outputFilename);
	void clearLog();
	void addLog(std::string _msg);
	void addLog(std::stringstream & _msg);


	std::string getLog();

	void writeLogToDisk(std::stringstream & _msg);
	void writeLogToDisk();
	void appendLogToDisk();
};



inline Log::Log()
{
	outputFilename = "untitled.log";
	logMsg = "";
}

inline Log::~Log()
{
	outputFilename = "";
	logMsg = "";
}


inline void Log::setOutputFilename(std::string _outputFilename)
{
	outputFilename = _outputFilename;
}


inline void Log::clearLog()
{
	logMsg = "";
}

inline void Log::addLog(std::string _msg)
{
	logMsg += _msg;
}


inline void Log::addLog(std::stringstream & _msg)
{
	addLog(_msg.str());
	_msg.str("");
}


inline void Log::writeLogToDisk(std::stringstream & _msg)
{
	addLog(_msg);
	writeLogToDisk();
}


inline void Log::writeLogToDisk()
{
	std::ofstream outputFile( outputFilename.c_str(), std::ios::out);
	outputFile << logMsg;
	outputFile.close();
}

inline void Log::appendLogToDisk()
{
	std::ofstream outputFile( outputFilename.c_str(), std::ios::out | std::ios::app);
	outputFile << logMsg;
	logMsg = "";
	outputFile.close();
}


inline std::string Log::getLog()
{
	return logMsg;
}



///////////////////////////////////////////////////////////////////////////////////
///////////// Simple Logging

inline void writeFile(std::string filename, std::string log)
{
	std::ofstream outputFile( filename.c_str(), std::ios::out);
	outputFile << log;
	outputFile.close();
}



inline void writeLog(std::string filename, std::string log)
{
	std::ofstream outputFile( (filename+ ".log").c_str(), std::ios::out);
	outputFile << log;
	outputFile.close();
}

inline void writeLogApp(std::string filename, std::string log)
{
	std::ofstream outputFile( (filename+ ".log").c_str(), std::ios::out | std::ios::app);
	outputFile << log;
	outputFile.close();
}



inline void writeLog(std::string filename, std::stringstream log)
{
	std::ofstream outputFile( (filename+ ".log").c_str(), std::ios::out);
	outputFile << log.str();
	outputFile.close();
}

inline void appendLog(std::string filename, std::stringstream & log)
{
	std::ofstream outputFile( (filename+ ".log").c_str(), std::ios::out | std::ios::app);
	outputFile << log.str();
	outputFile.close();

	log.str("");	// clears the log
}

#endif

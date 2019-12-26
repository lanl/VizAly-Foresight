/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/
#pragma once

#include <chrono>
#include <string>
#include <ctime>
#include <sstream>
#include <map>

class Timers
{
	std::map<
        std::string,
        std::chrono::time_point<std::chrono::system_clock> 
                > timers;
    std::map< 
        std::string, 
        std::chrono::duration<double> 
                > timers_duration;

  public:
	Timers();
	~Timers();

	int start(std::string timerName);
	int stop(std::string timerName);

	double getCurrentDuration(std::string timerName);	// time in seconds since timer started
	double getDuration(std::string timerName);		    // time in seconds

	static std::string getCurrentTime();	// get the current time
};

inline Timers::Timers() {}
inline Timers::~Timers() {}


inline int Timers::start(std::string timerName) 
{ 
	auto startTime = std::chrono::system_clock::now();
	if (timers.find(timerName) == timers.end())
		timers.insert( std::pair<std::string,std::chrono::time_point<std::chrono::system_clock>>(timerName,startTime) );	
	else
		timers[timerName] = startTime;

	return 1;
}


inline int Timers::stop(std::string timerName) 
{ 
	if (timers.find(timerName) != timers.end())
	{
		auto endTime = std::chrono::system_clock::now(); 
		auto elapsed_seconds = endTime - timers[timerName]; 
   		timers_duration.insert( std::pair<std::string,std::chrono::duration<double>>(timerName, elapsed_seconds) );
		   return 1;
	}
	else
	{
		std::cout << "Timer " << timerName << " does NOT exist!!!" << std::endl;
		return -1;
	}
}


inline double Timers::getDuration(std::string timerName) 
{ 
	if (timers_duration.find(timerName) != timers_duration.end())
	{

		return (timers_duration[timerName]).count(); 
	}
	else
	{
		std::cout << "Timer " << timerName << " does NOT exist!!!" << std::endl;
		return -1;
	}
}


inline double Timers::getCurrentDuration(std::string timerName)
{ 
	if (timers.find(timerName) != timers.end())
	{
		auto timeNow = std::chrono::system_clock::now();
		return (timeNow - timers[timerName]).count(); 
	}
	else
	{
		std::cout << "Timer " << timerName << " does NOT exist!!!" << std::endl;
		return -1;
	}
}


inline std::string Timers::getCurrentTime()
{
	time_t now = time(0);
	tm *ltm = localtime(&now);

	std::stringstream ss;
	ss << "_" << 1 + ltm->tm_mon << "_" << ltm->tm_mday << "__" << ltm->tm_hour << "_" << ltm->tm_min << "_" << ltm->tm_sec << "_" << std::endl;
	return ss.str();
}
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
	std::multimap<
        std::string,
        std::chrono::time_point<std::chrono::system_clock> 
                > timers;
    std::multimap< 
        std::string, 
        std::chrono::duration<double> 
                > timers_duration;

  public:
	Timers();
	~Timers();

	void start(std::string timerName);
	void stop(std::string timerName);

	double getCurrentDuration(timerName);	// time in seconds since timer started
	double getDuration(timerName);		    // time in seconds

	static std::string getCurrentTime();	// get the current time
};

inline Timers::Timers() {}
inline Timers::~Timers() {}


inline void Timers::start(std::string timerName) 
{ 
	auto startTime = std::chrono::system_clock::now(); 
	timers.insert( std::pair<std::string,std::chrono::time_point<std::chrono::system_clock>>(timerName,startTime) );
}


inline void Timers::stop(std::string timerName) 
{ 
	auto endTime = std::chrono::system_clock::now(); 
	auto elapsed_seconds = endTime - timers[timerName]; 
    timers_duration.insert( std::pair<std::string,std::chrono::duration<double>>(timerName, elapsed_seconds) );
}


inline double Timers::getDuration(std::string timerName) 
{ 
	return (timers_duration[timerName]).count(); 
}


inline double Timers::getCurrentDuration(std::string timerName)
{ 
	auto timeNow = std::chrono::system_clock::now();
	return (timeNow - timers[timerName]).count(); 
}


inline std::string Timers::getCurrentTime()
{
	time_t now = time(0);
	tm *ltm = localtime(&now);

	std::stringstream ss;
	ss << "_" << 1 + ltm->tm_mon << "_" << ltm->tm_mday << "__" << ltm->tm_hour << "_" << ltm->tm_min << "_" << ltm->tm_sec << "_" << std::endl;
	return ss.str();
}
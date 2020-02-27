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
		std::chrono::time_point<std::chrono::system_clock>
				> timers_duration;

  public:
  	enum TimeUnit{ sec=0, milli_sec=1, micro_sec=2, nano_sec=3 };

	Timers();
	~Timers();

	int start(std::string timerName);
	int stop(std::string timerName);

	double getDuration(std::string timerName, TimeUnit unit=TimeUnit::sec);	
};


inline Timers::Timers() {}
inline Timers::~Timers() {}


inline int Timers::start(std::string timerName) 
{ 
	auto startTime = std::chrono::system_clock::now();
	if (timers.find(timerName) == timers.end())
		timers.insert( std::pair<std::string, std::chrono::time_point<std::chrono::system_clock>>(timerName,startTime) );	
	else
		timers[timerName] = startTime;

	return 1;
}


inline int Timers::stop(std::string timerName) 
{ 
	if (timers.find(timerName) != timers.end())
	{
		auto endTime = std::chrono::system_clock::now(); 
		auto elapsed_time = endTime - timers[timerName];
		timers_duration.insert( std::pair<std::string,std::chrono::time_point<std::chrono::system_clock>>(timerName, elapsed_time) );
		return 1;
	}
	else
	{
		std::cout << "Timer " << timerName << " does NOT exist!!!" << std::endl;
		return -1;
	}
}


inline double Timers::getDuration(std::string timerName, TimeUnit unit) 
{ 
	if (timers_duration.find(timerName) != timers_duration.end())
	{
		if (unit == Timers::milli_sec)
			return ( std::chrono::duration_cast<std::chrono::milliseconds>(timers_duration[timerName]).count() );
		else if (unit == Timers::micro_sec)
			return ( std::chrono::duration_cast<std::chrono::microseconds>(timers_duration[timerName]).count() );
		else if (unit == Timers::nano_sec)
			return ( std::chrono::duration_cast<std::chrono::nanoseconds>(timers_duration[timerName]).count() ); 
		else
			return ( std::chrono::duration_cast<std::chrono::seconds>(timers_duration[timerName]).count() ); 
	}
	else
	{
		std::cout << "Timer " << timerName << " does NOT exist!!!" << std::endl;
		return -1;
	}
}



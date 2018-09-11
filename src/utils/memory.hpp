/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/

#ifndef _MEM_H_
#define _MEM_H_

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/sysinfo.h>
#include <unistd.h>
#endif // Linux

#if defined(WIN32)
#include <windows.h>
#include "psapi.h" //MSVC Performance query
#endif // Windows

#include <stdio.h>
#include <iostream>

class Memory
{
	unsigned long before_size, usage_size;	// program size
	unsigned long before_rss,  usage_rss;	// resident set size

	void GetMemorySize(unsigned long &size, unsigned long &rss);

  public:
	Memory();
	~Memory() {};

	void start();
	void stop();

	unsigned long getMemorySizeInB() { return usage_size; }
	double getMemorySizeInKB() { return usage_size / 1024.0; }
	double getMemorySizeInMB() { return usage_size / (1024.0 * 1024.0); }

	double getMemoryInUseInB();
	double getMemoryInUseInKB();
	double getMemoryInUseInMB();


	unsigned long getMemoryRSSInB() { return usage_rss; }
	double getMemoryRSSInKB() { return usage_rss / 1024.0; }
	double getMemoryRSSInMB() { return usage_rss / (1024.0 * 1024.0); }
};


inline Memory::Memory()
{
	before_size = usage_size = 0;
	before_rss  = usage_rss  = 0;
}


inline void Memory::start()
{
	GetMemorySize(before_size, before_rss);
}


inline void Memory::stop()
{
	unsigned long after_size, after_rss;
	GetMemorySize(after_size, after_rss);

	usage_size = after_size - before_size;
	usage_rss = after_rss - before_rss;
}


inline double Memory::getMemoryInUseInB()
{
	unsigned long after_size, after_rss;
	GetMemorySize(after_size, after_rss);

	return (after_size - before_size);
}

inline double Memory::getMemoryInUseInKB()
{
	unsigned long after_size, after_rss;
	GetMemorySize(after_size, after_rss);

	return (after_size - before_size) / (1024.0);
}


inline double Memory::getMemoryInUseInMB()
{
	unsigned long after_size, after_rss;
	GetMemorySize(after_size, after_rss);

	return (after_size - before_size) / (1024.0 * 1024.0);
}


#if defined(__unix__) || defined(__unix) || defined(unix)

// From VisIt avt/Pipeline/Pipeline/avtMemory.cpp
inline void Memory::GetMemorySize(unsigned long &size, unsigned long &rss)
{
	size = 0;
	rss  = 0;

	FILE *file = fopen("/proc/self/statm", "r");
	if (file == NULL)
		return;

	int count = fscanf(file, "%lu%lu", &size, &rss);
	if (count != 2)
	{
		fclose(file);
		return;
	}
	size *= (unsigned long)getpagesize();
	rss  *= (unsigned long)getpagesize();
	fclose(file);
}
#endif	// Linux

#if defined(WIN32)
inline void Memory::GetMemorySize(unsigned long &size, unsigned long &rss)
{
	//Virtual Memory by current process
	PROCESS_MEMORY_COUNTERS_EX pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
	SIZE_T virtualMemUsedByMe = pmc.PrivateUsage;

	//Physical Memory Used by Current Process
	SIZE_T physMemUsedByMe = pmc.WorkingSetSize;

	size = virtualMemUsedByMe;
	rss = physMemUsedByMe;
}
#endif // Windows

#endif	// _MEM_H_2

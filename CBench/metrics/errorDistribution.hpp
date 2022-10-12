/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/

#pragma once

#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include "metricInterface.hpp"


class ErrorDistribution : public MetricInterface
{
	int numRanks;
	int myRank;

public:
	ErrorDistribution();
	~ErrorDistribution();

	void init(MPI_Comm _comm);
	void execute(void *original, void *approx, size_t n, std::string dataType="float");
	void close() { }

};


inline ErrorDistribution::ErrorDistribution()
{
	myRank = 0;
	numRanks = 0;
	metricName = "error_distribution";
}


inline ErrorDistribution::~ErrorDistribution()
{

}


inline void ErrorDistribution::init(MPI_Comm _comm)
{
	comm = _comm;
	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);
}



inline void ErrorDistribution::ErrorDistribution(void *original, void *approx, size_t n, std::string dataType) {

	// Compute min max
	
	return;
}

#endif
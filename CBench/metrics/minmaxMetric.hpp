/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Jesus Pulido
================================================================================*/

#ifndef _MINMAX_METRIC_H_
#define _MINMAX_METRIC_H_

#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include "metricInterface.hpp"

class minmaxMetric : public MetricInterface
{
	int numRanks;
	int myRank;

public:
	minmaxMetric();
	~minmaxMetric();

	void init(MPI_Comm _comm);
	void execute(void *original, void *approx, size_t n);
	void close() { }

};

inline minmaxMetric::minmaxMetric()
{
	myRank = 0;
	numRanks = 0;
	metricName = "minmax";
}

inline minmaxMetric::~minmaxMetric()
{

}

inline void minmaxMetric::init(MPI_Comm _comm)
{
	comm = _comm;
	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);
}

inline void minmaxMetric::execute(void *original, void *approx, size_t n) {

	double local_max = -99999999999;
	double local_min = 99999999999;

	for (std::size_t i = 0; i < n; ++i)
	{
		if (static_cast<float *>(original)[i] > local_max)
			local_max = static_cast<float *>(original)[i];

		if (static_cast<float *>(original)[i] < local_min)
			local_min = static_cast<float *>(original)[i];

	}

	// Global max value
	double global_max = 0;
	MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm);

	// Global min value
	double global_min = 0;
	MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, comm);
    
	log << " local_minmax: " << local_min << " " << local_max << std::endl;
	// Currently only report Global minmax
	log << "-minmax: " << global_min << " " << global_max << std::endl;

	MPI_Barrier(comm);
	return;
}

#endif
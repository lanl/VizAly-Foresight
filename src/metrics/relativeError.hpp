/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/

#ifndef _RELATIVE_ERROR_H_
#define _RELATIVE_ERROR_H_

#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include "metricInterface.hpp"

class relativeError : public MetricInterface
{
	int numRanks;
	int myRank;

public:
	relativeError();
	~relativeError();

	void init(MPI_Comm _comm);
	void execute(void *original, void *approx, size_t n);
	void close() { }

};

inline relativeError::relativeError()
{
	myRank = 0;
	numRanks = 0;
	metricName = "relative_error";
}

inline relativeError::~relativeError()
{

}

inline void relativeError::init(MPI_Comm _comm)
{
	comm = _comm;
	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);
}


template <class T>
inline T relError(T original, T approx, double tolerance)
{
	double absolute_error = std::abs(original - approx);
	if (std::abs(original) < tolerance)
	{
		return absolute_error;
	}

	return absolute_error / std::abs(original);
}

inline void relativeError::execute(void *original, void *approx, size_t n) {
	std::vector<double> rel_err(n);

	for (std::size_t i = 0; i < n; ++i)
	{
		// Max set tolerence to 1
		double err = relError(static_cast<float *>(original)[i], static_cast<float *>(approx)[i], 1);
		rel_err.push_back(err);
	}
	double max_rel_err = *std::max_element(rel_err.begin(), rel_err.end());
	val = max_rel_err;

	double total_max_rel_err = 0;
	MPI_Reduce(&max_rel_err, &total_max_rel_err, 1, MPI_DOUBLE, MPI_MAX, 0, comm);// MPI_COMM_WORLD);
	total_val = total_max_rel_err;

	// Only reports the error per-rank for now
	log << " Max Rel Error: " << max_rel_err << std::endl;

	MPI_Barrier(comm);
	return;
}

#endif
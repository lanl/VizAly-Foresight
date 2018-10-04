/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/

#ifndef _ABSOLUTE_ERROR_H_
#define _ABSOLUTE_ERROR_H_

#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include "metricInterface.hpp"

class absoluteError : public MetricInterface
{
	int numRanks;
	int myRank;

public:
	absoluteError();
	~absoluteError();

	void init(MPI_Comm _comm);
	void execute(void *original, void *approx, size_t n);
	void close() { }

};


inline absoluteError::absoluteError()
{
	myRank = 0;
	numRanks = 0;
	metricName = "absolute_error";
}

inline absoluteError::~absoluteError()
{

}

inline void absoluteError::init(MPI_Comm _comm)
{
	comm = _comm;
	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);
}

template <class T>
inline T absError(T original, T approx)
{
	return std::abs(original - approx);
}

inline void absoluteError::execute(void *original, void *approx, size_t n) {
	std::vector<double> abs_err(n);

	for (std::size_t i = 0; i < n; ++i)
	{
		// Max set tolerence to 1
		double err = absError(static_cast<float *>(original)[i], static_cast<float *>(approx)[i]);
		abs_err.push_back(err);
	}
	double max_abs_err = *std::max_element(abs_err.begin(), abs_err.end());
	val = max_abs_err;

	double total_max_abs_err = 0;
	MPI_Allreduce(&max_abs_err, &total_max_abs_err, 1, MPI_DOUBLE, MPI_MAX, comm);// MPI_COMM_WORLD);
	total_val = total_max_abs_err;

	log << " Max Abs Error: " << total_max_abs_err << std::endl;

	MPI_Barrier(comm);
	return;
}

#endif
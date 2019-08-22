/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Jesus Pulido
================================================================================*/

#ifndef _PSNR_ERROR_H_
#define _PSNR_ERROR_H_

#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include "metricInterface.hpp"


class psnrError : public MetricInterface
{
	int numRanks;
	int myRank;

public:
	psnrError();
	~psnrError();

	void init(MPI_Comm _comm);
	void execute(void *original, void *approx, size_t n);
	void close() { }

};

inline psnrError::psnrError()
{
	myRank = 0;
	numRanks = 0;
	metricName = "psnr";
}

inline psnrError::~psnrError()
{

}

inline void psnrError::init(MPI_Comm _comm)
{
	comm = _comm;
	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);
}


/*template <class T>
inline T relError(T original, T approx, double tolerance)
{
	double absolute_error = std::abs(original - approx);
	if (std::abs(original) < tolerance)
	{
		return absolute_error;
	}

	return absolute_error / std::abs(original);
}*/

inline void psnrError::execute(void *original, void *approx, size_t n) {
	//std::vector<double> rel_err(n);

	double local_max = -999999999;
	double local_mse = 0;
	for (std::size_t i = 0; i < n; ++i)
	{
		if (static_cast<float *>(original)[i] > local_max)
			local_max = static_cast<float *>(original)[i];

		local_mse += (pow(static_cast<float *>(original)[i] - static_cast<float *>(approx)[i], (double)2.0));

		//double err;// = relError(static_cast<float *>(original)[i], static_cast<float *>(approx)[i], 1);
		//rel_err.push_back(err);
	}

	// Local Quantity
	double local_psnr = 10 * log10(pow(local_max, 2.0) / (local_mse / n));
	val = local_psnr;

	// Global max value
	double global_max = 0;
	MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm);

	// Global sum diff
	double global_mse = 0;
	MPI_Allreduce(&local_mse, &global_mse, 1, MPI_DOUBLE, MPI_SUM, comm);

	// Global number of values
	size_t global_n = 0;
	MPI_Allreduce(&n, &global_n, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

	global_mse /= global_n;
	double global_psnr = 10 * log10(pow(global_max, 2.0) / (global_mse));

	total_val = global_psnr;
	
	log << " local_psnr: " << local_psnr << std::endl;
	// Currently only report Global PSNR
	log << "-psnr: " << global_psnr << std::endl;

	MPI_Barrier(comm);
	return;
}

#endif
/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Jesus Pulido
================================================================================*/

#ifndef _MEAN_SQUARE_ERROR_H_
#define _MEAN_SQUARE_ERROR_H_

#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include "metricInterface.hpp"

class meansquareError : public MetricInterface
{
	int numRanks;
	int myRank;

public:
	meansquareError();
	~meansquareError();

	void init(MPI_Comm _comm);
	void execute(void *original, void *approx, size_t n, std::string dataType="float");
	void close() { }

};

inline meansquareError::meansquareError()
{
	myRank = 0;
	numRanks = 0;
	metricName = "mean_square_error";
}

inline meansquareError::~meansquareError()
{

}

inline void meansquareError::init(MPI_Comm _comm)
{
	comm = _comm;
	MPI_Comm_size(comm, &numRanks);
	MPI_Comm_rank(comm, &myRank);
}


inline void meansquareError::execute(void *original, void *approx, size_t n, std::string dataType) {
	//std::vector<double> mse(n);
	double mse = 0;

	for (std::size_t i = 0; i < n; ++i)
	{
		if (dataType == "float")
			mse += pow(( static_cast<float *>(original)[i]-static_cast<float *>(approx)[i]), 2.0);
		else if (dataType == "double")
			mse += pow(( static_cast<double *>(original)[i]-static_cast<double *>(approx)[i]), 2.0);
	}

	double local_mse = mse / n;
	val = local_mse;

	double total_mse = 0;
	size_t total_n = 0;

	MPI_Allreduce(&mse, &total_mse, 1, MPI_DOUBLE, MPI_SUM, comm);
	MPI_Allreduce(&n, &total_n, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
	total_val = total_mse/(double)total_n;
	
	debugLog << "-mse: " << total_val << std::endl;
	results["mse"] = mse;


	MPI_Barrier(comm);
	return;
}

#endif
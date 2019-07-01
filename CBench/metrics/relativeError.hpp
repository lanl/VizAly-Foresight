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

    double sum_rel_err = 0;
	for (std::size_t i = 0; i < n; ++i)
	{
		// Max set tolerence to 1
		double err = relError(static_cast<float *>(original)[i], static_cast<float *>(approx)[i], 1);
		rel_err.push_back(err);
        sum_rel_err += err;
	}
	double max_rel_err = *std::max_element(rel_err.begin(), rel_err.end());
	val = max_rel_err;

	double total_max_rel_err = 0;
	MPI_Allreduce(&max_rel_err, &total_max_rel_err, 1, MPI_DOUBLE, MPI_MAX, comm);// MPI_COMM_WORLD);
	total_val = total_max_rel_err;

	log << "-Max Rel Error: " << total_max_rel_err << std::endl;

    // Additional debug metrics, only in run_log
    // Global total sum of error
    double glob_sum_rel_err = 0;
    MPI_Allreduce(&sum_rel_err, &glob_sum_rel_err, 1, MPI_DOUBLE, MPI_SUM, comm);

    // Global number of values
    size_t global_n = 0;
    MPI_Allreduce(&n, &global_n, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

    // Compute mean
    double mean_rel_err = glob_sum_rel_err / global_n;
    log << " Total Rel Error: " << glob_sum_rel_err << std::endl;
    log << " Mean Rel Error: " << mean_rel_err << std::endl;

	MPI_Barrier(comm);


	auto found  = parameters.find("histogram");
	if ( found != parameters.end() )
	{
		// Compute histogram of values
		if (total_max_rel_err != 0)
		{
			std::vector<float>histogram;
			int numBins = 1024;
			std::vector<int> localHistogram(numBins,0);
			double binSize = total_max_rel_err / numBins;

			for (std::size_t i = 0; i < n; ++i)
			{
				// Max set tolerence to 1
				double err= relError(static_cast<float *>(original)[i], static_cast<float *>(approx)[i], 1);
				int binPos = err/binSize;

				if (binPos >= numBins)
					binPos = binPos-1;

				localHistogram[binPos]++;
			}

			histogram.resize(numBins);

			std::vector<int> globalHistogram(numBins,0);
			MPI_Allreduce(&localHistogram[0], &globalHistogram[0], numBins, MPI_INT, MPI_SUM, comm);

			for (std::size_t i=0; i<numBins; ++i)
				histogram[i] = ((float)globalHistogram[i])/global_n;


			// Output histogram as a python script file
			if (myRank == 0)
				additionalOutput = python_histogram(numBins,0.0, total_max_rel_err, histogram);
		}
	}


	return;
}

#endif
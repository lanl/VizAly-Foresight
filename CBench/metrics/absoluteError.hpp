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
#include <iostream>
#include <vector>
#include <math.h>

#include "metricInterface.hpp"

class absoluteError : public MetricInterface
{
	int numRanks;
	int myRank;

public:
	absoluteError();
	~absoluteError();

	void init(MPI_Comm _comm);
	void execute(void *original, void *approx, size_t n, std::string dataType="float");
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

inline void absoluteError::execute(void *original, void *approx, size_t n, std::string dataType) 
{
	std::vector<double> abs_err(n);

    double sum_abs_err = 0;
	for (std::size_t i = 0; i < n; ++i)
	{
		double err;
		if (dataType == "float")
			err = absError(static_cast<float *>(original)[i], static_cast<float *>(approx)[i]);
		else if (dataType == "double")
			err = absError(static_cast<double *>(original)[i], static_cast<double *>(approx)[i]);
			
		abs_err.push_back(err);
        sum_abs_err += err;
	}
	double max_abs_err = *std::max_element(abs_err.begin(), abs_err.end());
	val = max_abs_err;


	double total_max_abs_err = 0;
	MPI_Allreduce(&max_abs_err, &total_max_abs_err, 1, MPI_DOUBLE, MPI_MAX, comm);// MPI_COMM_WORLD);
	total_val = total_max_abs_err;

	debugLog << "-Local  Max Abs Error: " << max_abs_err << std::endl;
	debugLog << " Global Max Abs Error: " << total_max_abs_err << std::endl;
	results["absolute_error"] = total_max_abs_err;
	

    // Additional debug metrics, only in run_log
    // Global total sum of error
    double glob_sum_abs_err = 0;
    MPI_Allreduce(&sum_abs_err, &glob_sum_abs_err, 1, MPI_DOUBLE, MPI_SUM, comm);

    // Global number of values
    size_t global_n = 0;
    MPI_Allreduce(&n, &global_n, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

    // Compute mean
    double mean_abs_err = glob_sum_abs_err / global_n;
	double local_mean_abs_err = sum_abs_err/n;
    //debugLog << " Total Abs Error: " << glob_sum_abs_err << std::endl;
	debugLog << " Local Mean Abs Error : " << local_mean_abs_err << std::endl;
    debugLog << " Global Mean Abs Error: " << mean_abs_err << std::endl;
	results["mean_absolute_error"] = mean_abs_err;

	MPI_Barrier(comm);


	// Compute local standard deviation
	double localSumErrorSq = 0;
	for (std::size_t i = 0; i < n; ++i)
	{
		double err;
		if (dataType == "float")
			err = absError(static_cast<float *>(original)[i], static_cast<float *>(approx)[i]);
		else if (dataType == "double")
			err = absError(static_cast<double *>(original)[i], static_cast<double *>(approx)[i]);
			
        localSumErrorSq += (err - local_mean_abs_err)*(err - local_mean_abs_err);
	}
	debugLog << " Local STD Deviation : " << sqrt(localSumErrorSq/n) << std::endl;
	



	// Compute Global standard deviation
	double localGlobalSumErrorSq = 0;
	double globalSumErrorSq = 0;
	for (std::size_t i = 0; i < n; ++i)
	{
		double err;
		if (dataType == "float")
			err = absError(static_cast<float *>(original)[i], static_cast<float *>(approx)[i]);
		else if (dataType == "double")
			err = absError(static_cast<double *>(original)[i], static_cast<double *>(approx)[i]);
			
        localGlobalSumErrorSq += (err - mean_abs_err)*(err - mean_abs_err);
	}

    MPI_Allreduce(&localGlobalSumErrorSq, &globalSumErrorSq, 1, MPI_DOUBLE, MPI_SUM, comm);
	debugLog << " Gloabl STD Deviation : " << sqrt(globalSumErrorSq/global_n) << std::endl;
	results["std_absolute_error"] = sqrt(globalSumErrorSq/global_n);



	auto found  = parameters.find("histogram");
	if ( found != parameters.end() )
	{
		// Compute histogram of values
		if (total_max_abs_err != 0)
		{
			std::vector<float>histogram;
			size_t numBins = 1024;
			std::vector<size_t> localHistogram(numBins,0);
			double binSize = total_max_abs_err / numBins;

			for (std::size_t i = 0; i < n; ++i)
			{
				// Max set tolerence to 1
				double err = absError(static_cast<float *>(original)[i], static_cast<float *>(approx)[i]);

				int binPos = err/binSize;

				if (binPos >= numBins)
					binPos = binPos-1;

				localHistogram[binPos]++;
			}

			histogram.resize(numBins);

			std::vector<size_t> globalHistogram(numBins,0);
			MPI_Allreduce(&localHistogram[0], &globalHistogram[0], numBins, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

			for (std::size_t i=0; i<numBins; ++i)
				histogram[i] = ((float)globalHistogram[i])/(float)global_n;


			// Output histogram as a python script file
			if (myRank == 0)
				additionalOutput = python_histogram(numBins, 0.0, total_max_abs_err, histogram);
		}
	}
	

	return;
}

#endif
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
#include <limits>
#include "metricInterface.hpp"


class minmaxMetric : public MetricInterface
{
	int numRanks;
	int myRank;

public:
	minmaxMetric();
	~minmaxMetric();

	void init(MPI_Comm _comm);
	void execute(void *original, void *approx, size_t n, std::string dataType="float");
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


inline void minmaxMetric::execute(void *original, void *approx, size_t n, std::string dataType) {

	double local_max = -std::numeric_limits<double>::max();
	double local_min = std::numeric_limits<double>::max();

	for (std::size_t i = 0; i < n; ++i)
	{
		if (dataType == "float")
		{
			if (static_cast<float *>(original)[i] > local_max)
				local_max = static_cast<float *>(original)[i];

			if (static_cast<float *>(original)[i] < local_min)
				local_min = static_cast<float *>(original)[i];
		} 
		else if (dataType == "double")
		{
			if (static_cast<double *>(original)[i] > local_max)
				local_max = static_cast<double *>(original)[i];

			if (static_cast<double *>(original)[i] < local_min)
				local_min = static_cast<double *>(original)[i];
		}


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

	// Just report max for now
	val = local_max;
	total_val = global_max;

	auto found = parameters.find("histogram");
	if (found != parameters.end())
	{
		// Compute histogram of values
		if (global_max != 0)
		{
			std::vector<float>histogram;
			size_t numBins = 1024;
			std::vector<size_t> localHistogram(numBins, 0);
			double binSize = (global_max-global_min) / numBins;

			for (std::size_t i = 0; i < n; ++i)
			{
				// Retrieve the "approximated" value
				// Lossless comp: Original data
				// Lossy comp: Approx data
				double value = static_cast<float*>(approx)[i];

				int binPos = ( (global_max - global_min) * ((value - global_min)/ (global_max - global_min)) ) / binSize;

				if (binPos >= numBins)
					binPos = binPos - 1;

				localHistogram[binPos]++;
			}

			histogram.resize(numBins);

			std::vector<size_t> globalHistogram(numBins, 0);
			MPI_Allreduce(&localHistogram[0], &globalHistogram[0], numBins, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

			for (std::size_t i = 0; i < numBins; ++i)
				histogram[i] = (float) globalHistogram[i];

			// Output histogram as a python script file
			if (myRank == 0)
				additionalOutput = python_histogram(numBins, global_min, global_max, histogram);
				//additionalOutput = csv_histogram(numBins, global_max, histogram);
				
		}
	}

	return;
}

#endif
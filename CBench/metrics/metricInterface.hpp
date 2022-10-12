/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/

#ifndef _METRIC_INTERFACE_H_
#define _METRIC_INTERFACE_H_

#include <string>
#include <sstream>
#include <vector>
#include <mpi.h>

#include "log.hpp"


//
// synchronizes local hist with global histogram when using MPI, and returns it
inline std::vector<float> syncHistogram(int numBins, size_t numValues, std::vector<size_t> localHistogram, MPI_Comm comm)
{
    // Synchronize
    std::vector<size_t> globalHistogram(numBins, 0);
    MPI_Allreduce(&localHistogram[0], &globalHistogram[0], numBins, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

    // Get histogram
    std::vector<float>histogram;
    histogram.resize(numBins);

    for (std::size_t i=0; i<numBins; ++i)
        histogram[i] = ( (float)globalHistogram[i] )/( (float)numValues );
    
    return histogram;
}


class MetricInterface
{
  protected:
    double val;				// Local Quantity (MPI)
    double total_val;		// Global Quantity (MPI)

    std::string metricName;	// Internam metric name
    std::stringstream log;	// log file stream
    
	MPI_Comm comm;			// Global mpi comm handle

  public:
    std::unordered_map<std::string, double> results;
    std::unordered_map<std::string, std::string> parameters; // parameter inputs for metrics, i.e. enable histogram
    std::string additionalOutput;   // additional output stream, can be used for saving histogram (.py) file
                

	virtual void init(MPI_Comm _comm) = 0;
    virtual void execute(void *original, void *approx, size_t n, std::string dataType="float") = 0;
    virtual void close() = 0;

    std::string getMetricInfo();
    double getLocalValue(){ return val; }
    double getGlobalValue(){ return total_val; }
    std::string getMetricName(){ return metricName; }
    std::string getLog() { return log.str(); }
	void clearLog() { log.str(""); }
};



inline std::string MetricInterface::getMetricInfo()
{
    std::stringstream dataInfo;
    dataInfo << "\nMetric: " << metricName << std::endl;

    return dataInfo.str();
}

#endif

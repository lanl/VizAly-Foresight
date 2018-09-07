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

class MetricInterface
{
  protected:
    double val;
    double total_val;

    std::string metricName;
    std::stringstream log;

	MPI_Comm comm;
  public:
	virtual void init(MPI_Comm _comm) = 0;
    virtual void execute(void *original, void *approx, size_t n) = 0;
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

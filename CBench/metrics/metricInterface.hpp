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


inline std::string python_histogram(size_t numBins, float min_val, float max_val, std::vector<float> histogram)
{
    std::stringstream outputFileSS;
    outputFileSS << "import sys" << std::endl;
    outputFileSS << "import numpy as np" << std::endl;
	outputFileSS << "import matplotlib" << std::endl;
	outputFileSS << "matplotlib.use(\'agg\')" << std::endl;
	outputFileSS << "import matplotlib.pyplot as plt" << std::endl;

    outputFileSS << "y=[";
    std::size_t i;
    for (i=0; i<numBins-1; ++i) 
        outputFileSS << std::to_string(histogram[i]) << ", ";
    outputFileSS << std::to_string(histogram[i]) << "]" << std::endl;

	outputFileSS << "minVal=" << std::to_string(min_val) << std::endl;
    outputFileSS << "maxVal=" << std::to_string(max_val) << std::endl;
    outputFileSS << "plotName=sys.argv[0]" << std::endl;
    outputFileSS << "plotName = plotName.replace('.py','.png')" << std::endl;

    outputFileSS << "numVals = len(y)" << std::endl;
    outputFileSS << "x = np.linspace(minVal, maxVal, numVals+1)[1:]" << std::endl;
    outputFileSS << "plt.plot(x,y, linewidth=0.5)" << std::endl;
    outputFileSS << "plt.title(plotName)" << std::endl;
	outputFileSS << "plt.yscale(\"linear\") #log,linear,symlog,logit" << std::endl;
    outputFileSS << "plt.ylabel(\"Frequency\")" << std::endl;
    outputFileSS << "plt.xticks(rotation=90)" << std::endl;
    outputFileSS << "plt.tight_layout()" << std::endl;
    outputFileSS << "plt.savefig(plotName, dpi=300)" << std::endl;

    return outputFileSS.str();
}



class MetricInterface
{
  protected:
    double val;       // Local Quantity (MPI)
    double total_val; // Global Quantity (MPI)

    std::string metricName;
    std::stringstream log;
    
	MPI_Comm comm;

  public:
    std::unordered_map<std::string, std::string> parameters;
    std::string additionalOutput;   // if ever we need an additional output as for histograms
                

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

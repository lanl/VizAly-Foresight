/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/

#ifndef _METRICS_FACTORY_H_
#define _METRICS_FACTORY_H_

#include "metricIncludes.h"
#include "metricInterface.hpp"

class MetricsFactory
{
  public:
    static MetricInterface * createMetric(std::string metricName)
    {
        if (metricName == "absolute_error")
            return new absoluteError();
        else if (metricName == "relative_error")
            return new relativeError();
        else if (metricName == "mse")
            return new meansquareError();
        else
            return NULL;
    }
};

#endif
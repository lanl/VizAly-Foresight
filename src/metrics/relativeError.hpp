/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
 - Jesus Pulido
================================================================================*/

#pragma once
#include <cmath>

template <class T>
inline T relativeError(T original, T approx, double tolerance)
{
	double absolute_error = std::abs(original - approx);
	if (std::abs(original) < tolerance)
	{
		return absolute_error;
	}

	return absolute_error / std::abs(original);
}
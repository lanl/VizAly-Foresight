/*================================================================================
This software is open source software available under the BSD-3 license.

Copyright (c) 2017, Los Alamos National Security, LLC.
All rights reserved.

Authors:
 - Pascal Grosset
================================================================================*/

#pragma once
#include <cmath>

template <class T>
inline T absoluteError(T original, T approx)
{
	return std::abs(original - approx);
}
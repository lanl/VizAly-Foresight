#pragma once
#include <cmath>

template <class T>
inline T absoluteError(T original, T approx)
{
	return std::abs(original - approx);
}
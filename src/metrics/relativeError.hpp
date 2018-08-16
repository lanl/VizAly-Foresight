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
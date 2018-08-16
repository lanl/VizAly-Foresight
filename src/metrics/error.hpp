#ifndef ERROR_HPP
#define ERROR_HPP

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <tuple>

#include <mpi.h>




class error
{
public:
    error(const MPI_Comm &comm);

    template<class TYPE> double absolute(const TYPE *orig, const TYPE *pred, std::size_t size);
    template<class TYPE> std::tuple<double, double> relative(const TYPE *orig, const TYPE *pred, std::size_t size);
    template<class TYPE> double snr(const TYPE *orig, const TYPE *pred, std::size_t size);

private:
    const MPI_Comm &comm;

    template<class TYPE> double var(const TYPE *data, std::size_t size);
};


template<class TYPE> double error::absolute(const TYPE *orig, const TYPE *pred, std::size_t size)
{
    double diff = 0;
    for (std::size_t i = 0; i < size; ++i)
    {
        diff = std::max(diff, std::abs(static_cast<double>(pred[i]) - static_cast<double>(orig[i])));
    }

    double result;
    MPI_Allreduce(&diff, &result, 1, MPI_DOUBLE, MPI_MAX, this->comm);

    return result;
}


template<class TYPE> std::tuple<double, double> error::relative(const TYPE *orig, const TYPE *pred, std::size_t size)
{
    double rel_diff = 0;
    double abs_diff = 0;
    for (std::size_t i = 0; i < size; ++i)
    {
        double abs_err = std::abs(static_cast<double>(pred[i]) - static_cast<double>(orig[i]));
        if (orig[i] < 1)
        {
            abs_diff = std::max(abs_diff, abs_err);
        }
        else
        {
            rel_diff = std::max(rel_diff, abs_err / std::abs(static_cast<double>(orig[i])));
        }
    }

    double abs_result, rel_result;
    MPI_Allreduce(&abs_diff, &abs_result, 1, MPI_DOUBLE, MPI_MAX, this->comm);
    MPI_Allreduce(&rel_diff, &rel_result, 1, MPI_DOUBLE, MPI_MAX, this->comm);

    return std::make_tuple(abs_result, rel_result);
}


template<class TYPE> double error::snr(const TYPE *orig, const TYPE *pred, std::size_t size)
{
    std::vector<double> diff(size);
    for (std::size_t i = 0; i < size; ++i)
    {
        diff[i] = static_cast<double>(pred[i]) - static_cast<double>(orig[i]);
    }

    double result = 10.0 * std::log10(this->var(orig, size) / this->var(diff.data(), size));

    return result;
}


template<class TYPE> double error::var(const TYPE *data, std::size_t size)
{
    std::size_t tsize;
    MPI_Allreduce(&size, &tsize, 1, MPI_UINT64_T, MPI_SUM, this->comm);

    double temp = 0;
    for (std::size_t i = 0; i < size; ++i)
    {
        temp += static_cast<double>(data[i]);
    }

    double avg;
    MPI_Allreduce(&temp, &avg, 1, MPI_DOUBLE, MPI_SUM, this->comm);
    avg /= static_cast<double>(tsize);

    temp = 0;
    for (std::size_t i = 0; i < size; ++i)
    {
        temp += std::pow(data[i] - avg, 2);
    }

    double result;
    MPI_Allreduce(&temp, &result, 1, MPI_DOUBLE, MPI_SUM, this->comm);

    return result / static_cast<double>(tsize);
}


#endif // ERROR_HPP
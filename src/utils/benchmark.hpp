#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

#include <cstddef>
#include <vector>

#include <mpi.h>

class benchmark
{
public:

    benchmark(const MPI_Comm &comm);
    virtual ~benchmark() = default;

    void run();


protected:
    const MPI_Comm &comm;

    virtual void preprocess() = 0;
    virtual void init() = 0;
    virtual void executeComp() = 0;
    virtual void executeDecomp() = 0;
    virtual void cleanup() = 0;
    virtual void postprocess() = 0;
};

#endif // BENCHMARK_HPP
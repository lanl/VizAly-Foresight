#include <benchmark.hpp>

#include <algorithm>
#include <numeric>

#include <timer.hpp>




benchmark::benchmark(const MPI_Comm &comm)
    : comm(comm)
{

}


benchmark::timings benchmark::run()
{
    this->preprocess();

    //watch.trigger();
    this->init();
    //init_times[i] = watch.trigger();

    //watch.trigger();
    this->executeComp();
    //execute_times[i] = watch.trigger();

    this->executeDecomp();

    //watch.trigger();
    this->cleanup();
    //cleanup_times[i] = watch.trigger();

    this->postprocess();

    return;
}

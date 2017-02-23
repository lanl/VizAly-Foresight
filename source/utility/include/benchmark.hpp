#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

#include <cstddef>
#include <vector>

#include <mpi.h>




class benchmark
{
	public:
		struct results
		{
			double min_time, max_time, avg_time;
		};

		struct timings
		{
			benchmark::results init, execute, cleanup;
		};

		benchmark(const MPI_Comm &comm, std::size_t repetitions=1);
		virtual ~benchmark() = default;

		benchmark::timings run();

		std::size_t & repetitions();
		std::size_t repetitions() const;

	protected:
		const MPI_Comm &comm;
		std::size_t repetitions_m;

		virtual void preprocess() = 0;
		virtual void init() = 0;
		virtual void execute() = 0;
		virtual void cleanup() = 0;
		virtual void postprocess() = 0;
};


#endif // BENCHMARK_HPP

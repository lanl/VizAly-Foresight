#ifndef COMPRESSION_BENCHMARK_HPP
#define COMPRESSION_BENCHMARK_HPP

#include <vector>

extern "C"
{
#include <sz.h>
}

#include <data_field.hpp>
#include <benchmark.hpp>




class compression_benchmark : public benchmark
{
	public:
		std::vector<data_field> *data;
		std::vector<data_field> *cdata;

		sz_params params;

		compression_benchmark(const MPI_Comm &comm, std::size_t repetitions);
		~compression_benchmark();

	protected:
		void preprocess() override;
		void init() override;
		void execute() override;
		void cleanup() override;
		void postprocess() override;
};


#endif // COMPRESSION_BENCHMARK_HPP

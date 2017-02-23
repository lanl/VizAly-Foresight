#ifndef DECOMPRESSION_BENCHMARK_HPP
#define DECOMPRESSION_BENCHMARK_HPP

#include <vector>

#include <fpzip.h>

#include <data_field.hpp>
#include <benchmark.hpp>




class decompression_benchmark : public benchmark
{
	public:
		std::vector<data_field> *data;
		std::vector<data_field> *cdata;

		int precision;

		decompression_benchmark(const MPI_Comm &comm, std::size_t repetitions);
		~decompression_benchmark();

	private:
		std::vector<FPZ *> fpz;

	protected:
		void preprocess() override;
		void init() override;
		void execute() override;
		void cleanup() override;
		void postprocess() override;
};


#endif // DECOMPRESSION_BENCHMARK_HPP

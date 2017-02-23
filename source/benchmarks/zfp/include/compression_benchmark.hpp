#ifndef COMPRESSION_BENCHMARK_HPP
#define COMPRESSION_BENCHMARK_HPP

#include <vector>

#include <zfp.h>

#include <data_field.hpp>
#include <benchmark.hpp>




class compression_benchmark : public benchmark
{
	public:
		enum mode_t
		{
			ACCURACY,
			PRECISION,
			RATE
		};

		std::vector<data_field> *data;
		std::vector<data_field> *cdata;

		compression_benchmark::mode_t mode;
		std::size_t param;

		compression_benchmark(const MPI_Comm &comm, std::size_t repetitions);
		~compression_benchmark();

	private:
		std::vector<zfp_stream *> zfp;
		std::vector<zfp_field *> field;
		std::vector<bitstream *> stream;

	protected:
		void preprocess() override;
		void init() override;
		void execute() override;
		void cleanup() override;
		void postprocess() override;
};


#endif // COMPRESSION_BENCHMARK_HPP

#ifndef DECOMPRESSION_BENCHMARK_HPP
#define DECOMPRESSION_BENCHMARK_HPP

#include <vector>

#include <zfp.h>

#include <data_field.hpp>
#include <benchmark.hpp>




class decompression_benchmark : public benchmark
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

		decompression_benchmark::mode_t mode;
		std::size_t param;

		decompression_benchmark(const MPI_Comm &comm, std::size_t repetitions);
		~decompression_benchmark();

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


#endif // DECOMPRESSION_BENCHMARK_HPP

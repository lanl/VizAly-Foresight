#ifndef COMPRESSION_BENCHMARK_HPP
#define COMPRESSION_BENCHMARK_HPP

#include <cstdint>
#include <string>

#include <lzma.h>

#include <data_field.hpp>
#include <benchmark.hpp>




class compression_benchmark : public benchmark
{
	public:
		std::vector<data_field> *data;
		std::vector<data_field> *cdata;

		std::string preset;

		compression_benchmark(const MPI_Comm &comm, std::size_t repetitions);
		~compression_benchmark();

	private:
		std::uint32_t clevel;
		std::vector<lzma_stream> strm;
		std::uint8_t inbuf[BUFSIZ], outbuf[BUFSIZ];

	protected:
		void preprocess() override;
		void init() override;
		void execute() override;
		void cleanup() override;
		void postprocess() override;
};


#endif // COMPRESSION_BENCHMARK_HPP

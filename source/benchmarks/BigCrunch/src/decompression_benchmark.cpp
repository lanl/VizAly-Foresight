#include <decompression_benchmark.hpp>

#include <stdexcept>
#include <cstdlib>

extern "C"
{
#include <blosc.h>
}




decompression_benchmark::decompression_benchmark(const MPI_Comm &comm, std::size_t repetitions)
	: benchmark(comm, repetitions),
	data(nullptr),
	cdata(nullptr)
{
	this->settings = {{bigcrunch::config_t::ERROR,                                             -3},
                          {bigcrunch::config_t::TOLERANCE,                                          0},
                          {bigcrunch::config_t::BLOSC_NTHREADS,                                     1},
                          {bigcrunch::config_t::BLOSC_FILTER,      bigcrunch::blosc_filter_t::SHUFFLE},
                          {bigcrunch::config_t::BLOSC_COMPRESSOR, bigcrunch::blosc_compressor_t::ZSTD}};
}


decompression_benchmark::~decompression_benchmark()
{
	this->data = nullptr;
	this->cdata = nullptr;
}


void decompression_benchmark::preprocess()
{
	for(std::size_t i = 0; i < this->cdata->size(); ++i)
	{
		if(this->data->at(i).data != nullptr)
		{
			std::free(this->data->at(i).data);
		}
	}
}


void decompression_benchmark::init()
{
}


void decompression_benchmark::execute()
{
	bigcrunch::bigcrunch bc(this->settings);

	for(std::size_t i = 0; i < this->cdata->size(); ++i)
	{
		auto rdata_array = bc.decompress(static_cast<std::uint8_t *>(this->cdata->at(i).data), this->cdata->at(i).size);
		this->data->at(i).data = rdata_array.data<float>();
		this->data->at(i).size = rdata_array.size();
		this->data->at(i).type = data_field::data_type::FLOAT;
	}
}


void decompression_benchmark::cleanup()
{
}


void decompression_benchmark::postprocess()
{
}

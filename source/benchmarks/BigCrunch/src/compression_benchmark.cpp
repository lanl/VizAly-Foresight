#include <compression_benchmark.hpp>

#include <stdexcept>
#include <cstdlib>
#include <cstdint>
#include <algorithm>




compression_benchmark::compression_benchmark(const MPI_Comm &comm, std::size_t repetitions)
	: benchmark(comm, repetitions),
	data(nullptr),
	cdata(nullptr)
{
	this->settings = {{bigcrunch::config_t::ERROR,						 	  -3},
		         {bigcrunch::config_t::TOLERANCE,						   0},
		         {bigcrunch::config_t::BLOSC_NTHREADS,	 				 	   1},
		    	 {bigcrunch::config_t::BLOSC_FILTER,		 bigcrunch::blosc_filter_t::SHUFFLE},
		    	 {bigcrunch::config_t::BLOSC_COMPRESSOR,	 bigcrunch::blosc_compressor_t::ZSTD}};
}


compression_benchmark::~compression_benchmark()
{
	this->data = nullptr;
	this->cdata = nullptr;
}


void compression_benchmark::preprocess()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		if(this->cdata->at(i).data != nullptr)
		{
			std::free(this->cdata->at(i).data);
		}
	}
}


void compression_benchmark::init()
{
}


void compression_benchmark::execute()
{
	bigcrunch::bigcrunch bc(this->settings);

	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		std::uint8_t *cdata = nullptr;
		std::uint64_t csize = bc.compress(bigcrunch::darray(static_cast<float *>(this->data->at(i).data), this->data->at(i).size), &cdata);
		this->cdata->at(i).type = data_field::data_type::BYTE;
		this->cdata->at(i).size = csize;
		this->cdata->at(i).data = cdata;
	}
}


void compression_benchmark::cleanup()
{
}


void compression_benchmark::postprocess()
{
}

#include <compression_benchmark.hpp>

#include <stdexcept>
#include <cstdlib>
#include <algorithm>

extern "C"
{
#include <blosc.h>
}




compression_benchmark::compression_benchmark(const MPI_Comm &comm, std::size_t repetitions)
	: benchmark(comm, repetitions),
	data(nullptr),
	cdata(nullptr),
	nthreads(1),
	clevel(9),
	shuffle(1),
	cname("blosclz")
{
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
	blosc_init();
	blosc_set_nthreads(this->nthreads);
	int code = blosc_set_compressor(this->cname.c_str());
	if(code < 0)
	{
		throw std::invalid_argument("Error setting " + this->cname + " compressor. Does it really exist?");
	}
}


void compression_benchmark::execute()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		std::size_t typesize;
		switch(this->data->at(i).type)
		{
			case data_field::data_type::FLOAT:
				typesize = sizeof(float);
				break;
			case data_field::data_type::DOUBLE:
				typesize = sizeof(double);
				break;
			default:
				throw std::runtime_error("Invalid data type!");
		}

		std::size_t isize = typesize * this->data->at(i).size;
		std::size_t osize = isize + BLOSC_MAX_OVERHEAD;
		this->cdata->at(i).type = data_field::data_type::BYTE;
		this->cdata->at(i).size = isize;
		this->cdata->at(i).data = std::malloc(isize);

		this->cdata->at(i).size = blosc_compress(this->clevel, this->shuffle, typesize, isize, this->data->at(i).data, this->cdata->at(i).data, osize);
		if(this->cdata->at(i).size < 0)
		{
			throw std::runtime_error("Compression error. Error code: " + std::to_string(this->cdata->at(i).size));
		}

		if(this->cdata->at(i).size > 0)
		{
			this->cdata->at(i).data = std::realloc(this->cdata->at(i).data, this->cdata->at(i).size);
		}
	}
}


void compression_benchmark::cleanup()
{
	blosc_destroy();
}


void compression_benchmark::postprocess()
{
}

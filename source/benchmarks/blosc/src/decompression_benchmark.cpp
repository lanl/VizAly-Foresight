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
	cdata(nullptr),
	nthreads(1),
	clevel(9),
	shuffle(1),
	cname("blosclz")
{
}


decompression_benchmark::~decompression_benchmark()
{
	this->data = nullptr;
	this->cdata = nullptr;
}


void decompression_benchmark::preprocess()
{
}


void decompression_benchmark::init()
{
	blosc_init();
	blosc_set_nthreads(this->nthreads);
}


void decompression_benchmark::execute()
{
	for(std::size_t i = 0; i < this->cdata->size(); ++i)
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

		std::size_t osize = typesize * this->data->at(i).size;
		this->data->at(i).data = std::malloc(osize);
		this->data->at(i).size = blosc_decompress(this->cdata->at(i).data, this->data->at(i).data, osize);

		if(this->data->at(i).size < 0)
		{
			throw std::runtime_error("Decompression error. Error code: " + std::to_string(this->data->at(i).size));
		}
	}
}


void decompression_benchmark::cleanup()
{
	blosc_destroy();
}


void decompression_benchmark::postprocess()
{
}

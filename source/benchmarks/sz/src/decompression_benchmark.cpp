#include <decompression_benchmark.hpp>

#include <stdexcept>
#include <cstdint>
#include <algorithm>




decompression_benchmark::decompression_benchmark(const MPI_Comm &comm, std::size_t repetitions)
	: benchmark(comm, repetitions),
	data(nullptr),
	cdata(nullptr)
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
	SZ_Init_Params(&this->params);
}


void decompression_benchmark::execute()
{
	for(std::size_t i = 0; i < this->cdata->size(); ++i)
	{
		switch(this->data->at(i).type)
		{
			case data_field::data_type::FLOAT:
				this->data->at(i).data = static_cast<float *>(SZ_decompress(SZ_FLOAT, static_cast<std::uint8_t *>(this->cdata->at(i).data), this->cdata->at(i).size, 0, 0, 0, 0, this->data->at(i).size));
				break;
			case data_field::data_type::DOUBLE:
				this->data->at(i).data = static_cast<double *>(SZ_decompress(SZ_DOUBLE, static_cast<std::uint8_t *>(this->cdata->at(i).data), this->cdata->at(i).size, 0, 0, 0, 0, this->data->at(i).size));
				break;
			default:
				throw std::runtime_error("Invalid data type!");
		}
	}
}


void decompression_benchmark::cleanup()
{
	SZ_Finalize();
}


void decompression_benchmark::postprocess()
{
}

#include <compression_benchmark.hpp>

#include <stdexcept>
#include <algorithm>
#include <cstdint>




compression_benchmark::compression_benchmark(const MPI_Comm &comm, std::size_t repetitions)
	: benchmark(comm, repetitions),
	data(nullptr),
	cdata(nullptr)
{
}


compression_benchmark::~compression_benchmark()
{
	this->data = nullptr;
	this->cdata = nullptr;
}


void compression_benchmark::preprocess()
{
}


void compression_benchmark::init()
{
	SZ_Init_Params(&this->params);
}


void compression_benchmark::execute()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		int csize;
		switch(this->data->at(i).type)
		{
			case data_field::data_type::FLOAT:
				this->cdata->at(i).data = SZ_compress(SZ_FLOAT, this->data->at(i).data, &csize, 0, 0, 0, 0, this->data->at(i).size);
				break;
			case data_field::data_type::DOUBLE:
				this->cdata->at(i).data = SZ_compress(SZ_DOUBLE, this->data->at(i).data, &csize, 0, 0, 0, 0, this->data->at(i).size);
				break;
			default:
				throw std::runtime_error("Invalid data type!");
		}
		this->cdata->at(i).size = csize;
		this->cdata->at(i).type = data_field::data_type::BYTE;
	}
}


void compression_benchmark::cleanup()
{
	SZ_Finalize();
}


void compression_benchmark::postprocess()
{
}

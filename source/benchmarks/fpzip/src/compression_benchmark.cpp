#include <compression_benchmark.hpp>

#include <stdexcept>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <cstdlib>




compression_benchmark::compression_benchmark(const MPI_Comm &comm, std::size_t repetitions)
	: benchmark(comm, repetitions),
	data(nullptr),
	cdata(nullptr),
	precision(32)
{
}


compression_benchmark::~compression_benchmark()
{
	this->data = nullptr;
	this->cdata = nullptr;
}


void compression_benchmark::preprocess()
{
	this->fpz.clear();

	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		std::size_t limit;
		switch(this->data->at(i).type)
		{
			case data_field::data_type::FLOAT:
				limit = sizeof(float) * 8 - std::numeric_limits<float>::digits + 1;
				break;
			case data_field::data_type::DOUBLE:
				limit = sizeof(double) * 8 - std::numeric_limits<double>::digits + 1;
				break;
			default:
				throw std::runtime_error("Invalid data type!");
		}

		if(this->precision < limit + 1)
		{
			throw std::invalid_argument("Precision is too low for this datatype! Precision needs to be above " + std::to_string(limit) + ".");
		}
	}
}


void compression_benchmark::init()
{
}


void compression_benchmark::execute()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		this->fpz.push_back(fpzip_write_to_buffer(this->cdata->at(i).data, this->cdata->at(i).size));
		this->fpz[i]->prec = this->precision;
		this->fpz[i]->nx = this->data->at(i).size;
		this->fpz[i]->ny = 1;
		this->fpz[i]->nz = 1;
		this->fpz[i]->nf = 1;
		switch(this->data->at(i).type)
		{
			case data_field::data_type::FLOAT:
				this->fpz[i]->type = FPZIP_TYPE_FLOAT;
				this->cdata->at(i).size = 1024 + this->data->at(i).size * sizeof(float);
				break;
			case data_field::data_type::DOUBLE:
				this->fpz[i]->type = FPZIP_TYPE_DOUBLE;
				this->cdata->at(i).size = 1024 + this->data->at(i).size * sizeof(double);
			default:
				throw std::runtime_error("Invalid data type!");
		}

		this->cdata->at(i).type = data_field::data_type::BYTE;
		this->cdata->at(i).data = std::malloc(this->cdata->at(i).size);

		this->cdata->at(i).size = fpzip_write(this->fpz[i], this->data->at(i).data);
		this->cdata->at(i).data = std::realloc(this->cdata->at(i).data, this->cdata->at(i).size);
	}
}


void compression_benchmark::cleanup()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		fpzip_write_close(this->fpz[i]);
	}
}


void compression_benchmark::postprocess()
{
}

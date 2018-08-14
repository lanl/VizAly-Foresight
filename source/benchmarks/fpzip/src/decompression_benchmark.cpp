#include <decompression_benchmark.hpp>

#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <limits>




decompression_benchmark::decompression_benchmark(const MPI_Comm &comm, std::size_t repetitions)
	: benchmark(comm, repetitions),
	data(nullptr),
	cdata(nullptr),
	precision(32)
{
}


decompression_benchmark::~decompression_benchmark()
{
	this->data = nullptr;
	this->cdata = nullptr;
}


void decompression_benchmark::preprocess()
{
	this->fpz.clear();

	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		if(this->data->at(i).data != nullptr)
		{
			std::free(this->data->at(i).data);
		}

		std::size_t limit;
		switch(this->data->at(i).type)
		{
			case data_field::data_type::FLOAT:
				limit = sizeof(float) * 8 - std::numeric_limits<float>::digits + 1;
				break;
			case data_field::data_type::DOUBLE:
				limit = sizeof(double) * 8 - std::numeric_limits<double>::digits + 1;
			default:
				throw std::runtime_error("Invalid data type!");
		}

		if(this->precision < limit + 1)
		{
			throw std::invalid_argument("Precision is too low for this datatype! Precision needs to be above " + std::to_string(limit) + ".");
		}
	}
}


void decompression_benchmark::init()
{
}


void decompression_benchmark::execute()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		this->fpz.push_back(fpzip_read_from_buffer(static_cast<std::uint8_t *>(this->cdata->at(i).data)));
		this->fpz[i]->prec = this->precision;
		this->fpz[i]->nx = this->data->at(i).size;
		this->fpz[i]->ny = 1;
		this->fpz[i]->nz = 1;
		this->fpz[i]->nf = 1;

		switch(this->data->at(i).type)
		{
			case data_field::data_type::FLOAT:
				this->fpz[i]->type = FPZIP_TYPE_FLOAT;
				this->data->at(i).data = std::malloc(sizeof(float) * this->data->at(i).size);
				fpzip_read(this->fpz[i], static_cast<float *>(this->data->at(i).data));
				break;
			case data_field::data_type::DOUBLE:
				this->fpz[i]->type = FPZIP_TYPE_DOUBLE;
				this->data->at(i).data = std::malloc(sizeof(double) * this->data->at(i).size);
				fpzip_read(this->fpz[i], static_cast<double *>(this->data->at(i).data));
				break;
			default:
				throw std::runtime_error("Invalid data type!");
		}
	}
}


void decompression_benchmark::cleanup()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		fpzip_read_close(this->fpz[i]);
	}
}


void decompression_benchmark::postprocess()
{
}

#include <decompression_benchmark.hpp>

#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <cstdlib>
#include <cmath>




decompression_benchmark::decompression_benchmark(const MPI_Comm &comm, std::size_t repetitions)
	: benchmark(comm, repetitions),
	data(nullptr),
	cdata(nullptr),
	mode(decompression_benchmark::mode_t::ACCURACY),
	param(1.0 / std::pow(10, 6))
{
}


decompression_benchmark::~decompression_benchmark()
{
	this->data = nullptr;
	this->cdata = nullptr;
}


void decompression_benchmark::preprocess()
{
	this->zfp.clear();
	this->field.clear();
	this->stream.clear();
}


void decompression_benchmark::init()
{
}


void decompression_benchmark::execute()
{
	for(std::size_t i = 0; i < this->cdata->size(); ++i)
	{
		this->zfp.push_back(zfp_stream_open(nullptr));

		zfp_type type;
		switch(this->data->at(i).type)
		{
			case data_field::data_type::FLOAT:
				type = zfp_type_float;
				this->data->at(i).data = std::malloc(sizeof(float) * this->data->at(i).size);
				break;
			case data_field::data_type::DOUBLE:
				type = zfp_type_double;
				this->data->at(i).data = std::malloc(sizeof(double) * this->data->at(i).size);
				break;
			default:
				throw std::runtime_error("Invalid data type!");
		}

		this->stream.push_back(stream_open(this->cdata->at(i).data, this->cdata->at(i).size));
		zfp_stream_set_bit_stream(this->zfp[i], this->stream[i]);

		this->field.push_back(zfp_field_alloc());
		zfp_field_set_type(this->field[i], type);
		zfp_field_set_pointer(this->field[i], this->data->at(i).data);
		zfp_field_set_size_1d(this->field[i], this->data->at(i).size);

		switch(this->mode)
		{
			case decompression_benchmark::mode_t::ACCURACY:
				zfp_stream_set_accuracy(this->zfp[i], 1.0 / std::pow(10, this->param), type);
				break;
			case decompression_benchmark::mode_t::PRECISION:
				zfp_stream_set_precision(this->zfp[i], this->param, type);
				break;
			case decompression_benchmark::mode_t::RATE:
				zfp_stream_set_rate(this->zfp[i], this->param, type, 1, 0);
				break;
			default:
				throw std::invalid_argument("Unknown mode!");
		}

		zfp_stream_rewind(this->zfp[i]);

		zfp_decompress(this->zfp[i], this->field[i]);
	}
}


void decompression_benchmark::cleanup()
{
	for(std::size_t i = 0; i < this->cdata->size(); ++i)
	{
		zfp_field_free(this->field[i]);
		zfp_stream_close(this->zfp[i]);
		stream_close(this->stream[i]);
	}
}


void decompression_benchmark::postprocess()
{
}

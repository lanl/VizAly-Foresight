#include <compression_benchmark.hpp>

#include <stdexcept>
#include <algorithm>
#include <cstdint>
#include <cstdlib>

#include <GenericIOCompress.h>




compression_benchmark::compression_benchmark(const MPI_Comm &comm, std::size_t repetitions)
	: benchmark(comm, repetitions),
	data(nullptr),
	cdata(nullptr),
	num_buckets(1),
	blocksize(1024)
{
}


compression_benchmark::~compression_benchmark()
{
	this->data = nullptr;
	this->cdata = nullptr;
}


void compression_benchmark::preprocess()
{
	this->temp_cdata.clear();
	this->temp_cdata.resize(this->data->size());
}


void compression_benchmark::init()
{
}


void compression_benchmark::execute()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		switch(this->data->at(i).type)
		{
			case data_field::data_type::FLOAT:
				gio::sld_compressor<float> sldc_float;
				sldc_float.compressData_omp(static_cast<float *>(this->data->at(i).data), static_cast<float *>(this->data->at(i).data) + this->data->at(i).size, this->num_buckets, this->blocksize, this->temp_cdata[i]);
				break;
			case data_field::data_type::DOUBLE:
				gio::sld_compressor<double> sldc_double;
				sldc_double.compressData_omp(static_cast<double *>(this->data->at(i).data), static_cast<double *>(this->data->at(i).data) + this->data->at(i).size, this->num_buckets, this->blocksize, this->temp_cdata[i]);
				break;
			default:
				throw std::runtime_error("Invalid data type!");
		}
	}
}


void compression_benchmark::cleanup()
{
}


void compression_benchmark::postprocess()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		this->cdata->at(i).type = data_field::data_type::BYTE;
		this->cdata->at(i).size = this->temp_cdata[i].size();
		this->cdata->at(i).data = std::malloc(this->cdata->at(i).size);
		std::copy(this->temp_cdata[i].begin(), this->temp_cdata[i].end(), static_cast<std::uint8_t *>(this->cdata->at(i).data));
	}
}

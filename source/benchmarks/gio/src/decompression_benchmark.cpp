#include <decompression_benchmark.hpp>

#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <cstdlib>

#include <GenericIOCompress.h>




decompression_benchmark::decompression_benchmark(const MPI_Comm &comm, std::size_t repetitions)
	: benchmark(comm, repetitions),
	data(nullptr),
	cdata(nullptr),
	num_buckets(1),
	blocksize(1024)
{
}


decompression_benchmark::~decompression_benchmark()
{
	this->data = nullptr;
	this->cdata = nullptr;
}


void decompression_benchmark::preprocess()
{
	this->temp_cdata.clear();
	this->temp_cdata.resize(this->cdata->size());

	for(std::size_t i = 0; i < this->cdata->size(); ++i)
	{
		this->temp_cdata[i].resize(this->cdata->at(i).size);
		std::copy(static_cast<std::uint8_t *>(this->cdata->at(i).data), static_cast<std::uint8_t *>(this->cdata->at(i).data) + this->cdata->at(i).size, this->temp_cdata[i].begin());
	}
}


void decompression_benchmark::init()
{
}


void decompression_benchmark::execute()
{
	for(std::size_t i = 0; i < this->cdata->size(); ++i)
	{
		switch(this->data->at(i).type)
		{
			case data_field::data_type::FLOAT:
				this->data->at(i).data = std::malloc(sizeof(float) * this->data->at(i).size);
				gio::sld_compressor<float> sldc_float;
				sldc_float.decompressData_omp(static_cast<float *>(this->data->at(i).data), static_cast<float *>(this->data->at(i).data) + this->data->at(i).size, this->num_buckets, this->blocksize, this->temp_cdata[i], 0);
				break;
			case data_field::data_type::DOUBLE:
				this->data->at(i).data = std::malloc(sizeof(double) * this->data->at(i).size);
				gio::sld_compressor<double> sldc_double;
				sldc_double.decompressData_omp(static_cast<double *>(this->data->at(i).data), static_cast<double *>(this->data->at(i).data) + this->data->at(i).size, this->num_buckets, this->blocksize, this->temp_cdata[i], 0);
				break;
			default:
				throw std::runtime_error("Invalid data type!");
		}
	}
}


void decompression_benchmark::cleanup()
{
}


void decompression_benchmark::postprocess()
{
}

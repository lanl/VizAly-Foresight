#include <compression_benchmark.hpp>

#include <stdexcept>
#include <algorithm>
#include <cstdlib>




compression_benchmark::compression_benchmark(const MPI_Comm &comm, std::size_t repetitions)
	: benchmark(comm, repetitions),
	data(nullptr),
	cdata(nullptr),
	preset("9e")
{
}


compression_benchmark::~compression_benchmark()
{
	this->data = nullptr;
	this->cdata = nullptr;
}


void compression_benchmark::preprocess()
{
	this->clevel = this->preset[0] - '0';
	if(this->preset.size() > 1)
	{
		if(this->preset[1] != 'e' || this->preset.size() > 2)
		{
			throw std::invalid_argument("Given preset is invalid!");
		}

		this->clevel |= LZMA_PRESET_EXTREME;
	}
}


void compression_benchmark::init()
{
	this->strm = LZMA_STREAM_INIT;
	lzma_ret ret = lzma_easy_encoder(&this->strm, this->clevel, LZMA_CHECK_CRC64);
	switch(ret)
	{
		case LZMA_OK:
			break;
		case LZMA_MEM_ERROR:
			throw std::runtime_error("Memory allocation failed!");
			break;
		case LZMA_OPTIONS_ERROR:
			throw std::runtime_error("Specified preset is not supported!");
			break;
		case LZMA_UNSUPPORTED_CHECK:
			throw std::runtime_error("Specified integrity check is not supported!");
			break;
		default:
			throw std::runtime_error("Unknown error, possibly a bug!");
	}
}


void compression_benchmark::execute()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		lzma_action action = LZMA_RUN;

		this->strm.next_in = nullptr;
		this->strm.avail_in = 0;
		this->strm.next_out = this->outbuf;
		this->strm.avail_out = sizeof(this->outbuf);

		std::size_t it = 0;
		bool running = true;

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
		this->cdata->at(i).size = 0;
		this->cdata->at(i).type = data_field::data_type::BYTE;
		this->cdata->at(i).data = std::malloc(typesize * this->data->at(i).size);

		while(running)
		{
			if(this->strm.avail_in == 0 && it < this->data->at(i).size)
			{
				this->strm.next_in = this->inbuf;
				std::size_t bufsize = 0;
				while(bufsize + typesize < BUFSIZ && it < this->data->at(i).size)
				{
					for(int i = 0; i < typesize; ++i)
					{
						this->inbuf[bufsize + i] = static_cast<std::uint8_t *>(this->data->at(i).data)[it];
					}
					bufsize += typesize;
					++it;
				}
				this->strm.avail_in = bufsize;

				if(it == this->data->at(i).size)
				{
					action = LZMA_FINISH;
				}
			}

			lzma_ret ret = lzma_code(&this->strm, action);
			if(this->strm.avail_out == 0 || ret == LZMA_STREAM_END)
			{
				std::size_t write_size = sizeof(this->outbuf) - this->strm.avail_out;
				std::copy(this->outbuf, this->outbuf + write_size, static_cast<std::uint8_t *>(this->cdata->at(i).data) + this->cdata->at(i).size);
				this->cdata->at(i).size += write_size;
				this->strm.next_out = this->outbuf;
				this->strm.avail_out = sizeof(this->outbuf);
			}

			switch(ret)
			{
				case LZMA_OK:
					break;
				case LZMA_STREAM_END:
					running = false;
					break;
				case LZMA_MEM_ERROR:
					throw std::runtime_error("Memory allocation failed!");
					break;
				case LZMA_DATA_ERROR:
					throw std::runtime_error("File size limits exceeded!");
					break;
				default:
					throw std::runtime_error("Unknown error, possibly a bug!");
			}
		}

		this->cdata->at(i).data = std::realloc(this->cdata->at(i).data, this->cdata->at(i).size);
	}
}


void compression_benchmark::cleanup()
{
	lzma_end(&this->strm);
}


void compression_benchmark::postprocess()
{
}

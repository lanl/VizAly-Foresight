#include <decompression_benchmark.hpp>

#include <stdexcept>
#include <algorithm>
#include <cstdlib>




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
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		if(this->data->at(i).data != nullptr)
		{
			std::free(this->data->at(i).data);
		}
	}
}


void decompression_benchmark::init()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		this->strm.push_back(LZMA_STREAM_INIT);
		lzma_ret ret = lzma_stream_decoder(&this->strm[i], UINT64_MAX, LZMA_CONCATENATED);
		switch(ret)
		{
			case LZMA_OK:
				break;
			case LZMA_MEM_ERROR:
				throw std::runtime_error("Memory allocation failed!");
				break;
			case LZMA_OPTIONS_ERROR:
				throw std::runtime_error("Unsupported decompressor flags!");
				break;
			default:
				throw std::runtime_error("Unknown error, possibly a buf!");
		}
	}
}


void decompression_benchmark::execute()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		lzma_action action = LZMA_RUN;

		this->strm[i].next_in = nullptr;
		this->strm[i].avail_in = 0;
		this->strm[i].next_out = this->outbuf;
		this->strm[i].avail_out = sizeof(this->outbuf);

		std::size_t it = 0;
		std::size_t pos = 0;
		bool running = true;

		switch(this->data->at(i).type)
		{
			case data_field::data_type::FLOAT:
				this->data->at(i).data = std::malloc(sizeof(float) * this->data->at(i).size);
				break;
			case data_field::data_type::DOUBLE:
				this->data->at(i).data = std::malloc(sizeof(float) * this->data->at(i).size);
				break;
			default:
				throw std::runtime_error("Invalid data type!");
		}

		while(running)
		{
			if(this->strm[i].avail_in == 0 && it < this->cdata->at(i).size)
			{
				this->strm[i].next_in = this->inbuf;
				std::size_t bufsize = 0;
				while(bufsize + 1 < BUFSIZ && it < this->cdata->at(i).size)
				{
					this->inbuf[bufsize] = static_cast<std::uint8_t *>(this->cdata->at(i).data)[it];
					++bufsize;
					++it;
				}
				this->strm[i].avail_in = bufsize;

				if(it == this->cdata->at(i).size)
				{
					action = LZMA_FINISH;
				}
			}

			lzma_ret ret = lzma_code(&this->strm[i], action);
			if(this->strm[i].avail_out == 0 || ret == LZMA_STREAM_END)
			{
				std::size_t write_size = sizeof(this->outbuf) - this->strm[i].avail_out;
				std::copy(this->outbuf, this->outbuf + write_size, static_cast<std::uint8_t *>(this->data->at(i).data) + pos);
				pos += write_size;
				this->strm[i].next_out = this->outbuf;
				this->strm[i].avail_out = sizeof(outbuf);
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
				case LZMA_FORMAT_ERROR:
					throw std::runtime_error("The input is not in the .xz format!");
					break;
				case LZMA_OPTIONS_ERROR:
					throw std::runtime_error("Unsupported compression options!");
					break;
				case LZMA_DATA_ERROR:
					throw std::runtime_error("Compressed file is corrupt!");
					break;
				case LZMA_BUF_ERROR:
					throw std::runtime_error("Compressed file is truncated or otherwise corrupt!");
					break;
				default:
					throw std::runtime_error("Unknown error, possibly a bug!");
			}
		}
	}
}


void decompression_benchmark::cleanup()
{
	for(std::size_t i = 0; i < this->data->size(); ++i)
	{
		lzma_end(&this->strm[i]);
	}
}


void decompression_benchmark::postprocess()
{
}

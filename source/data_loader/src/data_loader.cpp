#include <data_loader.hpp>

#include <stdexcept>

#include <hacc_data.hpp>
#include <mpas_data.hpp>




data_handler *  data_loader::create(const MPI_Comm &comm, data_loader::format data_format)
{
	switch(data_format)
	{
		case data_loader::format::HACC:
			return new hacc_data(comm);
		case data_loader::format::MPAS:
			return new mpas_data(comm);
		default:
			throw std::invalid_argument("Unknown data format!");
	}
}

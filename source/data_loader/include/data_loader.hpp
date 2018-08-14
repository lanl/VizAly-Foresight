#ifndef DATA_LOADER_HPP
#define DATA_LOADER_HPP

#include <string>

#include <mpi.h>

#include <data_handler.hpp>
#include <hacc_data.hpp>
#include <mpas_data.hpp>




class data_loader
{
	public:
		enum format
		{
			HACC,
			MPAS
		};

		static data_handler * create(MPI_Comm &comm, data_loader::format data_format);
};


#endif // DATA_LOADER_HPP

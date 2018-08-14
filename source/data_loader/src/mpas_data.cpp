#include <mpas_data.hpp>

#include <iostream>
#include <stdexcept>




mpas_data::mpas_data(const MPI_Comm &comm)
	: comm(comm)
{
}


void mpas_data::read(const std::string &in_name, std::vector<data_field> &data)
{
	data.clear();
	data.resize(2);

	int rank, size;
	MPI_Comm_rank(this->comm, &rank);
	MPI_Comm_size(this->comm, &size);

	data[0].type = data_field::data_type::DOUBLE;
	data[0].size = size;

	data[1].type = data_field::data_type::DOUBLE;
	data[1].size = size;
}


void mpas_data::write(const std::string &in_name, const std::string &out_name, const std::vector<data_field> &data)
{
	int rank;
	MPI_Comm_rank(this->comm, &rank);
}

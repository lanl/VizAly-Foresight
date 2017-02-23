#include <mpas_data.hpp>




mpas_data::mpas_data(const MPI_Comm &comm)
	: comm(comm)
{
}


void mpas_data::read(const std::string &in_name, std::vector<data_field> &data)
{
}


void mpas_data::write(const std::string &in_name, const std::string &out_name, const std::vector<data_field> &data)
{
}

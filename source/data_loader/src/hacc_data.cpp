#include <hacc_data.hpp>

#include <cstdlib>
#include <cstdint>

#include <GenericIO.h>




hacc_data::hacc_data(const MPI_Comm &comm)
	: comm(comm)
{
}


void hacc_data::read(const std::string &in_name, std::vector<data_field> &data)
{
	data.clear();
	std::vector<std::vector<float>> temp;
	std::vector<std::string> channel = {"x", "y", "z", "vx", "vy", "vz"};
	temp.resize(channel.size());

	int rank;
	MPI_Comm_rank(this->comm, &rank);

	gio::GenericIO handler(this->comm, in_name);
	handler.openAndReadHeader();
	std::size_t num_particles = handler.readNumElems(rank);
	for(int i = 0; i < temp.size(); ++i)
	{
		temp[i].resize(num_particles + handler.requestedExtraSpace() / sizeof(float));
		handler.addVariable(channel[i], temp[i], true);
	}
	handler.readData(rank);

	data.resize(channel.size());
	for(int i = 0; i < temp.size(); ++i)
	{
		data[i].type = data_field::data_type::FLOAT;
		data[i].size = temp[i].size();
		data[i].data = std::malloc(sizeof(float) * data[i].size);
		std::copy(temp[i].begin(), temp[i].end(), static_cast<float *>(data[i].data));
	}
}


void hacc_data::write(const std::string &in_name, const std::string &out_name, const std::vector<data_field> &data)
{
	//---- Data containers ----
	std::vector<float> x, y, z, vx, vy, vz, phi;
	std::vector<std::int64_t> id;
	std::vector<std::uint16_t> mask;

	//---- Get MPI information ----
	int rank;
	MPI_Comm_rank(this->comm, &rank);

	//---- Header information ----
	//int num_ranks;
	//int dims[3];
	double phys_origin[3], phys_scale[3];
	std::vector<gio::GenericIO::VariableInfo> vinfo;
	//std::vector<rankData> rank_meta;

	//---- Open GenericIO file and read header ----
	gio::GenericIO in_handler(this->comm, in_name);
	in_handler.openAndReadHeader();

	//---- Extract GenericIO header information ----
	std::size_t num_particles = in_handler.readNumElems(rank);
	//in_handler.readDims(dims);
	in_handler.readPhysOrigin(phys_origin);
	in_handler.readPhysScale(phys_scale);
	//num_ranks = in_handler.readNRanks();
	in_handler.getVariableInfo(vinfo);

	/*rank_meta.resize(num_ranks);
	for(std::size_t i = 0; i < rank_meta.size(); ++i)
	{
		rank_meta[i].numPoints = in_handler.readNumElems(i);
		in_handler.readCoords(rank_meta[i].coords, i);
	}*/	

	//---- Add variables to be read to the gio file handler ----
	x.resize(num_particles + in_handler.requestedExtraSpace() / sizeof(float));
	in_handler.addVariable("x", x, true);
	y.resize(num_particles + in_handler.requestedExtraSpace() / sizeof(float));
	in_handler.addVariable("y", y, true);
	z.resize(num_particles + in_handler.requestedExtraSpace() / sizeof(float));
	in_handler.addVariable("z", z, true);
	vx.resize(num_particles + in_handler.requestedExtraSpace() / sizeof(float));
	in_handler.addVariable("vx", vx, true);
	vy.resize(num_particles + in_handler.requestedExtraSpace() / sizeof(float));
	in_handler.addVariable("vy", vy, true);
	vz.resize(num_particles + in_handler.requestedExtraSpace() / sizeof(float));
	in_handler.addVariable("vz", vz, true);
	phi.resize(num_particles + in_handler.requestedExtraSpace() / sizeof(float));
	in_handler.addVariable("phi", phi, true);
	id.resize(num_particles + in_handler.requestedExtraSpace() / sizeof(std::int64_t));
	in_handler.addVariable("id", id, true);
	mask.resize(num_particles + in_handler.requestedExtraSpace() / sizeof(std::uint16_t));
	in_handler.addVariable("mask", mask, true);

	//---- Read data ----
	in_handler.readData(rank);

	//---- Replace particle data ----
	x.clear();
	x.assign(static_cast<float *>(data[0].data), static_cast<float *>(data[0].data) + data[0].size);
	y.clear();
	y.assign(static_cast<float *>(data[1].data), static_cast<float *>(data[1].data) + data[1].size);
	z.clear();
	z.assign(static_cast<float *>(data[2].data), static_cast<float *>(data[2].data) + data[2].size);
	vx.clear();
	vx.assign(static_cast<float *>(data[3].data), static_cast<float *>(data[3].data) + data[3].size);
	vy.clear();
	vy.assign(static_cast<float *>(data[4].data), static_cast<float *>(data[4].data) + data[4].size);
	vz.clear();
	vz.assign(static_cast<float *>(data[5].data), static_cast<float *>(data[5].data) + data[5].size);

	//---- Open GenericIO output file ----
	gio::GenericIO out_handler(this->comm, out_name);

	//---- Input GenericIO header information ----
	for(std::size_t i = 0; i < 3; ++i)
	{
		out_handler.setPhysOrigin(phys_origin[i], i);
		out_handler.setPhysScale(phys_scale[i], i);
	}

	out_handler.setNumElems(num_particles);

	//---- Add variables to be written to the gio file handler ----
	out_handler.addVariable(vinfo[0], x.data());
	out_handler.addVariable(vinfo[1], y.data());
	out_handler.addVariable(vinfo[2], z.data());
	out_handler.addVariable(vinfo[3], vx.data());
	out_handler.addVariable(vinfo[4], vy.data());
	out_handler.addVariable(vinfo[5], vz.data());
	out_handler.addVariable(vinfo[6], phi.data());
	out_handler.addVariable(vinfo[7], id.data());
	out_handler.addVariable(vinfo[8], mask.data());

	//---- Write data ----
	out_handler.write();
}

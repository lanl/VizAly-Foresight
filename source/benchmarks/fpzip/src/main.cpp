#include <sstream>
#include <iomanip>

#include <data_loader.hpp>
#include <compression_benchmark.hpp>
#include <decompression_benchmark.hpp>
#include <fio.hpp>
#include <error.hpp>




std::size_t gather_size(const MPI_Comm &comm, const data_field &field)
{
	std::size_t rank_size;
	switch(field.type)
	{
		case data_field::data_type::FLOAT:
			rank_size = sizeof(float) * field.size;
			break;
		case data_field::data_type::DOUBLE:
			rank_size = sizeof(double) * field.size;
			break;
		case data_field::data_type::BYTE:
			rank_size = field.size;
			break;
		default:
			throw std::runtime_error("Invalid data type!");
	}

	std::size_t total_size = 0;
	MPI_Reduce(&rank_size, &total_size, 1, MPI_UINT64_T, MPI_SUM, 0, comm);

	return total_size;
}


int main(int argc, char **argv)
{
	//---- Gather command line arguments ----
	int minargs = 7;
	if(argc - 1 < minargs)
	{
		throw std::runtime_error("Wrong number of command line arguments! Minimum " + std::to_string(minargs) + " were expected, " + std::to_string(argc - 1) + " were given.");
	}

	std::string in_name = std::string(argv[1]);
	std::string out_name = std::string(argv[2]);
	std::string result_name = std::string(argv[3]);
	std::string format = std::string(argv[4]);
	std::size_t repetitions = std::stoul(std::string(argv[5]));
	int precision = std::stoi(std::string(argv[6]));
	std::size_t nthreads = std::stoul(std::string(argv[7]));

	//---- Init MPI ----
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	//---- Load data ----
	data_handler *handler = nullptr;

	if(format.compare("HACC") == 0)
	{
		handler = data_loader::create(comm, data_loader::format::HACC);
	}

	if(format.compare("MPAS") == 0)
	{
		handler = data_loader::create(comm, data_loader::format::MPAS);
	}

	std::vector<data_field> data;
	handler->read(in_name, data);

	std::vector<data_field> cdata, udata;
	cdata.resize(data.size());
	udata.resize(data.size());
	for(std::size_t i = 0; i < data.size(); ++i)
	{
		udata[i].size = data[i].size;
		udata[i].type = data[i].type;
	}

	//---- Run compression benchmark ----
	compression_benchmark cbench(comm, repetitions);
	cbench.data = &data;
	cbench.cdata = &cdata;
	cbench.precision = precision;
	benchmark::timings cresults = cbench.run();

	//---- Run decompression benchmark ----
	decompression_benchmark dbench(comm, repetitions);
	dbench.data = &udata;
	dbench.cdata = &cdata;
	dbench.precision = precision;
	benchmark::timings dresults = dbench.run();

	//---- Save reconstructed data ----
	//if(nthreads == 16)
	//{
	//	handler->write(in_name, out_name, udata);
	//}

	//---- Gather total original and compressed data sizes ----
	std::vector<std::size_t> original_size(data.size() + 1, 0);
	std::vector<std::size_t> compressed_size(data.size() + 1, 0);
	for(std::size_t i = 0; i < data.size(); ++i)
	{
		original_size[i] = gather_size(comm, data[i]);
		compressed_size[i] = gather_size(comm, cdata[i]);
		original_size.back() += original_size[i];
		compressed_size.back() += compressed_size[i];
	}

	//---- Compute error metrics ----
	std::vector<double> abs_err, rel_err, snr;
	error err(comm);
	for(std::size_t i = 0; i < data.size(); ++i)
	{
		std::tuple<double, double> error;
		switch(data[i].type)
		{
			case data_field::data_type::FLOAT:
				error = err.relative(static_cast<float *>(data[i].data), static_cast<float *>(udata[i].data), data[i].size);
				abs_err.push_back(std::get<0>(error));
				rel_err.push_back(std::get<1>(error));
				snr.push_back(err.snr(static_cast<float *>(data[i].data), static_cast<float *>(udata[i].data), data[i].size));
				break;
			case data_field::data_type::DOUBLE:
				error = err.relative(static_cast<double *>(data[i].data), static_cast<double *>(udata[i].data), data[i].size);
				abs_err.push_back(std::get<0>(error));
				rel_err.push_back(std::get<1>(error));
				snr.push_back(err.snr(static_cast<double *>(data[i].data), static_cast<double *>(udata[i].data), data[i].size));
				break;
		}
	}

	//---- Save results ----
	if(rank == 0)
	{
		std::vector<double> ratios(data.size() + 1, 0);
		for(std::size_t i = 0; i < ratios.size(); ++i)
		{
			ratios[i] = original_size[i] / static_cast<double>(compressed_size[i]);
		}

		std::stringstream ss;
		ss << in_name << ";" << format << ";" << repetitions << ";" << size << ";" << precision << ";" << nthreads << ";";
		ss << original_size.back() << ";" << compressed_size.back() << ";";
		ss << std::setprecision(15) << ratios.back() << ";";

		for(std::size_t i = 0; i < data.size(); ++i)
		{
			ss << ratios[i] << ";";
		}

		for(std::size_t i = 0; i < data.size(); ++i)
		{
			ss << abs_err[i] << ";";
		}

		for(std::size_t i = 0; i < data.size(); ++i)
		{
			ss << rel_err[i] << ";";
		}

		for(std::size_t i = 0; i < data.size(); ++i)
		{
			ss << snr[i] << ";";
		}

		double gb = original_size.back() / 1000.0 / 1000.0 / 1000.0;

		ss << gb / cresults.init.min_time << ";" << gb / cresults.init.max_time << ";" << gb / cresults.init.avg_time << ";";
		ss << gb / cresults.execute.min_time << ";" << gb / cresults.execute.max_time << ";" << gb / cresults.execute.avg_time << ";";
		ss << gb / cresults.cleanup.min_time << ";" << gb / cresults.cleanup.max_time << ";" << gb / cresults.cleanup.avg_time << ";";

		ss << gb / dresults.init.min_time << ";" << gb / dresults.init.max_time << ";" << gb / dresults.init.avg_time << ";";
		ss << gb / dresults.execute.min_time << ";" << gb / dresults.execute.max_time << ";" << gb / dresults.execute.avg_time << ";";
		ss << gb / dresults.cleanup.min_time << ";" << gb / dresults.cleanup.max_time << ";" << gb / dresults.cleanup.avg_time;

		fio file(result_name, std::ios::out | std::ios::app);
		file.writeline(ss.str());
	}

	//---- Cleanup ----
	data.clear();
	cdata.clear();
	udata.clear();
	delete handler;
	MPI_Finalize();

	return 0;
}

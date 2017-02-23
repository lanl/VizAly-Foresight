#include <sstream>
#include <iomanip>

#include <data_loader.hpp>
#include <compression_benchmark.hpp>
#include <decompression_benchmark.hpp>
#include <fio.hpp>



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
	int minargs = 9;
	if(argc - 1 < minargs)
	{
		throw std::runtime_error("Wrong number of command line arguments! Minimum " + std::to_string(minargs) + " were expected, " + std::to_string(argc - 1) + " were given.");
	}

	std::string in_name = std::string(argv[1]);
	std::string out_name = std::string(argv[2]);
	std::string result_name = std::string(argv[3]);
	std::string format = std::string(argv[4]);
	std::size_t repetitions = std::stoul(std::string(argv[5]));
	int nthreads = std::stoi(std::string(argv[6]));
	int clevel = std::stoi(std::string(argv[7]));
	int shuffle = std::stoi(std::string(argv[8]));
	std::string cname = std::string(argv[9]);

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
	cbench.nthreads = nthreads;
	cbench.clevel = clevel;
	cbench.shuffle = shuffle;
	cbench.cname = cname;
	benchmark::timings cresults = cbench.run();

	//---- Run decompression benchmark ----
	decompression_benchmark dbench(comm, repetitions);
	dbench.data = &udata;
	dbench.cdata = &cdata;
	dbench.nthreads = nthreads;
	dbench.clevel = clevel;
	dbench.shuffle = shuffle;
	dbench.cname = cname;
	benchmark::timings dresults = dbench.run();

	//---- Save reconstructed data ----
	//handler->write(in_name, out_name, udata);

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

	//---- Save results ----
	if(rank == 0)
	{
		std::vector<double> ratios(data.size() + 1, 0);
		for(std::size_t i = 0; i < ratios.size(); ++i)
		{
			ratios[i] = original_size[i] / static_cast<double>(compressed_size[i]);
		}

		std::stringstream ss;
		ss << in_name << ";" << format << ";" << repetitions << ";"  << size << ";" << cname << ";" << clevel << ";" << shuffle << ";" << nthreads << ";";
		ss << original_size.back() << ";" << compressed_size.back() << ";";
		ss << std::setprecision(15) << ratios.back() << ";";

		for(std::size_t i = 0; i < data.size(); ++i)
		{
			ss << ratios[i] << ";";
		}

		ss << cresults.init.min_time << ";" << cresults.init.max_time << ";" << cresults.init.avg_time << ";";
		ss << cresults.execute.min_time << ";" << cresults.execute.max_time << ";" << cresults.execute.avg_time << ";";
		ss << cresults.cleanup.min_time << ";" << cresults.cleanup.max_time << ";" << cresults.cleanup.avg_time << ";";

		ss << dresults.init.min_time << ";" << dresults.init.max_time << ";" << dresults.init.avg_time << ";";
		ss << dresults.execute.min_time << ";" << dresults.execute.max_time << ";" << dresults.execute.avg_time << ";";
		ss << dresults.cleanup.min_time << ";" << dresults.cleanup.max_time << ";" << dresults.cleanup.avg_time << ";";

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

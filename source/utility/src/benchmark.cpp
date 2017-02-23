#include <benchmark.hpp>

#include <algorithm>
#include <numeric>

#include <stopwatch.hpp>




benchmark::benchmark(const MPI_Comm &comm, std::size_t repetitions)
	: comm(comm),
	repetitions_m(repetitions)
{
	if(repetitions < 1)
	{
		this->repetitions_m = 1;
	}
}


benchmark::timings benchmark::run()
{
	stopwatch watch(this->comm);

	std::vector<double> init_times(this->repetitions_m, 0);
	std::vector<double> execute_times(this->repetitions_m, 0);
	std::vector<double> cleanup_times(this->repetitions_m, 0);

	//for(std::size_t i = 0; i < this->repetitions_m; ++i)
	for(std::size_t i = 0; i < 1; ++i)
	{
		this->preprocess();

		//watch.trigger();
		this->init();
		//init_times[i] = watch.trigger();
		
		//watch.trigger();
		this->execute();
		//execute_times[i] = watch.trigger();

		//watch.trigger();
		this->cleanup();
		//cleanup_times[i] = watch.trigger();

		this->postprocess();
	}

	auto init_minmax = std::minmax_element(init_times.begin(), init_times.end());
	auto execute_minmax = std::minmax_element(execute_times.begin(), execute_times.end());
	auto cleanup_minmax = std::minmax_element(cleanup_times.begin(), cleanup_times.end());
	/*double avg = 0;
	for(auto &time : times)
	{
		avg += time;
	}
	avg /= static_cast<double>(times.size());*/

	benchmark::timings results = {
		.init = {
			.min_time = *(init_minmax.first),
			.max_time = *(init_minmax.second),
			.avg_time = std::accumulate(init_times.begin(), init_times.end(), 0.0) / static_cast<double>(init_times.size())
		},
		.execute = {
			.min_time = *(execute_minmax.first),
			.max_time = *(execute_minmax.second),
			.avg_time = std::accumulate(execute_times.begin(), execute_times.end(), 0.0) / static_cast<double>(execute_times.size())
		},
		.cleanup = {
			.min_time = *(cleanup_minmax.first),
			.max_time = *(cleanup_minmax.second),
			.avg_time = std::accumulate(cleanup_times.begin(), cleanup_times.end(), 0.0) / static_cast<double>(cleanup_times.size())
		}
	};

	/*benchmark::results result = {
		.min_time=*(minmax.first),
		.max_time=*(minmax.second),
		//.avg_time=avg
		.avg_time=std::accumulate(times.begin(), times.end(), 0.0) / static_cast<double>(times.size())
	};*/

	return results;
}


std::size_t & benchmark::repetitions()
{
	return this->repetitions_m;
}


std::size_t benchmark::repetitions() const
{
	return this->repetitions_m;
}

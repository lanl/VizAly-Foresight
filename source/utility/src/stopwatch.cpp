#include <stopwatch.hpp>




stopwatch::stopwatch(const MPI_Comm &comm)
	: comm(comm),
	last(MPI_Wtime())
{
}


double stopwatch::trigger()
{
	MPI_Barrier(this->comm);
	double current = MPI_Wtime();
	double time = current - this->last;
	this->last = current;

	return time;
}

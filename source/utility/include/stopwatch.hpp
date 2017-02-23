#ifndef STOPWATCH_HPP
#define STOPWATCH_HPP

#include <mpi.h>




class stopwatch
{
	public:
		stopwatch(const MPI_Comm &comm);

		double trigger();

	private:
		double last;
		const MPI_Comm &comm;
};


#endif // STOPWATCH_HPP

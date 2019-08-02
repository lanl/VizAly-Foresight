/* --------------------------------------------------------------------------
 * This software is open source software available under the BSD-3 license.
 *
 * Copyright (c) 2017, Los Alamos National Security, LLC.
 * All rights reserved.
 *
 * Author:
 * - Hoby Rakotoarivelo
 */
/* -------------------------------------------------------------------------- */
#include "analyzer.hpp"
/* -------------------------------------------------------------------------- */
int main(int argc, char* argv[]){

  int my_rank = 0;
  int nb_ranks = 0;
  int threading = 0;
  MPI_Comm comm = MPI_COMM_WORLD;
  
  // init MPI
  MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &threading);
  MPI_Comm_size(comm, &nb_ranks);
  MPI_Comm_rank(comm, &my_rank);

  // basic input checks
  if (not tools::valid(argc, argv, my_rank, nb_ranks)) {
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  // everything is OK at this point
  Analyzer analyzer(argv[1], my_rank, nb_ranks, comm);
 
  // run the analyzer 
  analyzer.run();

  MPI_Finalize();
  return EXIT_SUCCESS;
}
/* -------------------------------------------------------------------------- */

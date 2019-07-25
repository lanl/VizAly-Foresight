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
#include "halo_entropy_analyzer.hpp"

/* -------------------------------------------------------------------------- */
void printUsage() {
  std::fprintf(stderr, "Usage: mpirun -n <int> analyzer [options]\n");
  std::fprintf(stderr, "options:\n");
  std::fprintf(stderr, "-h, --help   show this help message and exit\n");
  std::fprintf(stderr, "-i, --input  input json parameter file\n");
  std::fflush(stderr);
}
/* -------------------------------------------------------------------------- */
bool isValidInput(int argc, char* argv[], int my_rank, int nb_ranks) {

  if (argc < 2) {
    if (my_rank == 0)
      printUsage();
    return false;
  } else if(not isPowerOfTwo(nb_ranks)) {
    std::cerr << "Please run with power of 2 ranks" << std::endl;
    return false;
  }

  // check input json file
  std::string path(argv[1]);
  std::ifstream file(path);
  nlohmann::json json;

  if (not file.good()) {
    if (my_rank == 0)
      std::cerr << "Error while opening parameter file: "<< argv[1] << std::endl;
    file.close();
    return false;
  }

  try {
    // pass file to json parser
    file >> json;

    // retrieve input file
    if (json["input"]["filetype"] != "HACC") {
      if (my_rank == 0)
        std::cerr << "Only HACC data is supported" << std::endl;
      file.close();
      return false;
    }

    file.close();
    return true;
  } catch(nlohmann::json::parse_error& e) {
    if (my_rank == 0)
      std::cerr << "Invalid JSON file " << path << "\n" << e.what() << std::endl;
    file.close();
    return false;
  }
}

/* -------------------------------------------------------------------------- */
int main(int argc, char* argv[]){

  int my_rank = 0;
  int nb_ranks = 0;
  int thread_support = 0;
  

  // init MPI
  MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &thread_support);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // basic input checks
  if (not isValidInput(argc, argv, my_rank, nb_ranks)) {
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  // everything is OK at this point
  HaloEntropy analyzer(argv[1], my_rank, nb_ranks, MPI_COMM_WORLD);
 
  // run the analyzer 
  analyzer.run();

  MPI_Finalize();
  return EXIT_SUCCESS;
}
/* -------------------------------------------------------------------------- */

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
#pragma once
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <mpi.h>

// helper functions
#include "json.hpp"
#include "timer.hpp"
#include "log.hpp"
#include "utils.hpp"
#include "strConvert.hpp"
/* -------------------------------------------------------------------------- */
namespace tools {
/* -------------------------------------------------------------------------- */
inline bool valid(int argc, char* argv[], int my_rank, int nb_ranks) {

  if (argc < 2) {
    if (my_rank == 0)
      std::cerr << "Usage: mpirun -n <int> ./analyzer [input-json]" << std::endl;
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
inline std::string base(std::string path) {
  char sep = '/';
#ifdef _WIN32
  sep = '\\';
#endif

  size_t i = path.rfind(sep, path.length());
  if (i != std::string::npos) {
    return std::string(path.substr(i+1, path.length() - i));
  }
  return std::string("");
}
/* -------------------------------------------------------------------------- */

} // namespace tools
/* -------------------------------------------------------------------------- */
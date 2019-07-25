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
#ifndef _HALO_ENTROPY_ANALYZER_H_
#define _HALO_ENTROPY_ANALYZER_H_
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

// data loader
#include "dataLoaderInterface.hpp"
#include "HACCDataLoader.hpp"

// EDIT: process one parameter each time
/* -------------------------------------------------------------------------- */
class HaloEntropy {

public: 
  HaloEntropy() = default; 
  HaloEntropy(HaloEntropy const&) = delete; 
  HaloEntropy(HaloEntropy&&) noexcept = delete; 
  HaloEntropy(const char* in_path, int in_rank, int in_nb_ranks, MPI_Comm in_comm);
  ~HaloEntropy() = default;
 
  bool run();  // ok
  bool computeFrequencies(std::string scalar, void* original, void* approx, size_t n);
  double computeShannonEntropy(std::string scalar, size_t n);
  bool generateHistogram(std::string path) const;

private:
  // IO
  std::string json_path = "";                                        // JSON parameter file path
  std::unique_ptr<DataLoaderInterface> ioMgr;                        // HACC data loader
  std::stringstream debug_log;

  // bins partition
  size_t const num_bins = 1024;
  // per halo attribute data
  std::vector<std::string> attributes;                               // [x,y,x,vx,vy,vx]
  std::unordered_map<std::string, std::vector<double>> frequency;
  std::unordered_map<std::string, std::string> type;                 // datatype of each attribute [float|double]
  std::unordered_map<std::string, int> size;                         // size of each attribute dataset (in bytes)
  std::unordered_map<std::string, int> nb_elems;                     // number of dataset elements for each attribute

  // MPI
  int my_rank  = 0;                                                   // current MPI rank
  int nb_ranks = 0;                                                  // total number of ranks
  MPI_Comm comm = MPI_COMM_WORLD;                                    // MPI communicator
};

/* -------------------------------------------------------------------------- */
inline HaloEntropy::HaloEntropy(const char* in_path, int in_rank,
                                int in_nb_ranks, MPI_Comm in_comm)
  : json_path(in_path), my_rank(in_rank),
    nb_ranks(in_nb_ranks), comm(in_comm),
    ioMgr(std::make_unique<HACCDataLoader>())
{}

/* -------------------------------------------------------------------------- */
inline bool HaloEntropy::computeFrequencies(std::string scalar, void* original, void* approx, size_t n) {

  // step 1. determine lower and upper bounds on data
  double global_max = 0;
  double global_min = 0;
  double local_max = std::numeric_limits<float>::min();
  double local_min = std::numeric_limits<float>::max();

  for (size_t i = 0; i < n; ++i) {
    if (static_cast<float*>(original)[i] > local_max)
      local_max = static_cast<float*>(original)[i];

    if (static_cast<float*>(original)[i] < local_min)
      local_min = static_cast<float*>(original)[i];
  }

  MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, comm);

  debug_log << " local_minmax: " << local_min << " " << local_max << std::endl;
  debug_log << " global_minmax: " << global_min << " " << global_max << std::endl;
  MPI_Barrier(comm);

  // Compute histogram of values
  if (global_max != 0) {

    debug_log << "num_bins: " << num_bins << std::endl;

    size_t local_histogram[num_bins];
    size_t global_histogram[num_bins];

    std::fill(local_histogram, local_histogram + num_bins, 0);
    std::fill(global_histogram, global_histogram + num_bins, 0);

    double range = global_max - global_min;
    double capacity = range / num_bins;

    for (size_t i = 0; i < n; ++i) {

      double value = static_cast<float*>(approx)[i];
      double relative_value = (value - global_min) / range;
      int index = (range * relative_value) / capacity;

      /*
      debug_log << "value: " << value
                << ", relative_value: " << relative_value
                << ", global_max: " << global_max << std::endl;*/
      if (index >= num_bins) {
        index--;
      }

      local_histogram[index]++;
    }

    MPI_Allreduce(local_histogram, global_histogram, num_bins, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

    // fill frequency eventually
    frequency[scalar].clear();
    frequency[scalar].resize(num_bins);

    for (size_t i = 0; i < num_bins; ++i) {
      frequency[scalar][i] = static_cast<double>(global_histogram[i]) / n;
      debug_log << "frequency[" << i <<"]: " << frequency[scalar][i] << std::endl;
    }
    return true;
  }

  return false;
}

/* -------------------------------------------------------------------------- */
inline double HaloEntropy::computeShannonEntropy(std::string scalar, size_t n) {

  auto const& distribution = frequency[scalar];

  double entropy = 0.;
  for (auto&& p_i : distribution) {
    entropy += (p_i * std::log2(p_i));
  }
  entropy *= -1;

  debug_log << "Shannon entropy for "<< scalar << ", using "
            << num_bins << " bins is " << entropy
            << ", accuracy: " << num_bins / static_cast<double>(n)
            << std::endl;

  return entropy;
}

/* -------------------------------------------------------------------------- */
inline bool HaloEntropy::run() {
  // basic checks 
  assert(nb_ranks > 0); 
  assert(ioMgr != nullptr);

  nlohmann::json json;
  std::ifstream file(json_path);

  debug_log.clear();
  debug_log.str("");

  try {
    // pass file to json parser
    file >> json;
  } catch(nlohmann::json::parse_error& e) {
    if (my_rank == 0) 
      std::cerr << "Invalid JSON file " << json_path << "\n" << e.what() << std::endl;
    return false;
  }

  // retrieve input file
  std::string filetype = json["input"]["filetype"];
  std::string input_hacc = json["input"]["filename"];
  if (filetype != "HACC") {
    if (my_rank == 0) 
      std::cerr << "Only HACC data is supported" << std::endl;
    return false;
  }

  // setup logs
  std::string basename = json["cbench"]["output"]["log-file"]; 
  std::string output = basename +"_rank_"+ std::to_string(my_rank) +".log";

  for (auto&& scalar : json["input"]["scalars"]) {
    attributes.push_back(scalar);
    // todo: retrieve metadata here
  }


  // set the IO manager
  if (json["input"].find("datainfo") != json["input"].end()) {
    auto const& datainfo = json["input"]["datainfo"];
    for (auto it = datainfo.begin(); it != datainfo.end(); ++it) {
      ioMgr->loaderParams[it.key()] = strConvert::toStr(it.value());
    }    
  }

  ioMgr->init(input_hacc, comm); 
  ioMgr->setSave(false); 
  MPI_Barrier(comm);

  if (my_rank == 0)
    printf("Found %d attributes\n", attributes.size());

  // loop over all scalars and load each of them
  for (auto&& scalar: attributes) {

    // for debug:
    if (scalar != "x") continue;

    debug_log << "\nLoading and running " << scalar << std::endl;
    if (ioMgr->loadData(scalar)) {
      // save infos for debug
      debug_log << ioMgr->getDataInfo();
      debug_log << ioMgr->getLog();
      // TODO: compute Shannon entropy distribution for this attribute
      if (scalar == "x") {
        auto const n = ioMgr->getNumElements();
        computeFrequencies(scalar, ioMgr->data, ioMgr->data, n);
        computeShannonEntropy(scalar, n);
      }

    } else {
      if (my_rank == 0)
        std::cerr << "Error while loading " << scalar << ", exiting now" << std::endl;
      return false;
    }
    //appendLog(output_log, log.str());
    // wait for other ranks to complete
    MPI_Barrier(comm);
  }

  // print some metadata on loaded HACC file

  // output logs
  std::ofstream logfile(output, std::ios::out);
  logfile << debug_log.str();
  debug_log.str("");
  logfile.close();

  std::cout << "Logs generated in "<< output << std::endl;  
  MPI_Barrier(comm);

  // everything was fine at this point  
  return true;
}

/* -------------------------------------------------------------------------- */


#endif

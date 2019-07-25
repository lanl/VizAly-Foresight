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
  bool computeFrequencies(std::string scalar, void* original, void* approx);
  double computeShannonEntropy(std::string scalar);
  bool generateHistogram(std::string path) const;
  void exportResults(std::string path = ".");
  std::string extractBase(std::string path);

private:
  // IO
  std::string json_path = "";
  std::string input_hacc = "";
  std::unique_ptr<DataLoaderInterface> ioMgr;
  std::stringstream debug_log;

  // bins partition
  size_t num_bins = 1024;
  // per halo attribute data
  std::vector<std::string> attributes;                               // [x,y,x,vx,vy,vx]
  std::unordered_map<std::string, std::vector<double>> frequency;
  std::unordered_map<std::string, size_t> size;                         // size of each attribute dataset (in bytes)

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
inline bool HaloEntropy::computeFrequencies(std::string scalar, void* original, void* approx) {

  double const& n = size[scalar];
  assert(n > 0.);

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

      if (index >= num_bins) { index--; }

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
inline double HaloEntropy::computeShannonEntropy(std::string scalar) {

  double entropy = 0.;
  double const epsilon = std::numeric_limits<double>::epsilon();

  for (auto&& value : frequency[scalar]) {
    auto const& p_i = std::max(value, epsilon);
    entropy += (p_i * std::log2(p_i));
  }
  entropy *= -1;

  debug_log << "shannon entropy for '"<< scalar << "' = "<< entropy
           << " using " << num_bins << " bins." << std::endl;

  return entropy;
}

/* -------------------------------------------------------------------------- */
inline std::string HaloEntropy::extractBase(std::string path) {
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
inline void HaloEntropy::exportResults(std::string root_path) {

  if (my_rank != 0)
    return;

  std::string const num_bins_str = std::to_string(num_bins);

  // dump data first
  for (auto&& scalar : attributes) {
    std::string path = root_path + "_" + scalar + "_" + num_bins_str +".dat";
    std::ofstream file(path, std::ios::out|std::ios::trunc);
    assert(file.is_open());
    assert(file.good());

    file << "# scalar: " << scalar << std::endl;
    file << "# num_bins: " << num_bins_str << std::endl;
    for (auto&& value : frequency[scalar]) {
      file << value << std::endl;
    }

    file.close();
  }

  // generate gnuplot script then
  std::string path = root_path + "_" + num_bins_str +".gnu";
  std::string output_eps = root_path +"_"+ num_bins_str +".eps";

  std::cout << path << std::endl;
  std::cout << output_eps << std::endl;

  std::ofstream file(path, std::ios::out);
  assert(file.is_open());
  assert(file.good());

  file << "#" << std::endl;
  file << "# HACC data distribution gnuplot script" << std::endl;
  file << "#" << std::endl;

  file << "reset" << std::endl;
  file << "set terminal postscript eps enhanced color 14 size 21cm,14cm" << std::endl;
  file << "set output '"<< output_eps <<"'" << std::endl;
  file << std::endl;
  file << "set multiplot layout 2,3 title 'HACC particle data - ";
  file << num_bins << " bins - file: " << extractBase(input_hacc) << "'" << std::endl;

  for (auto&& scalar : attributes) {
    // get bounds
    double v_max = 0.;
    for (auto&& value : frequency[scalar]) {
      if (value > v_max) { v_max = value; }
    }

    double const entropy = computeShannonEntropy(scalar);

    std::string data_file = root_path + "_" + scalar + "_" + num_bins_str +".dat";

    file << std::endl;
    file << "# ----------------------------------------------------" << std::endl;
    file << "set title '"<< scalar << ", entropy=" << entropy <<"'"<< std::endl;
    file << "set size ratio 1" << std::endl;
    file << "set xlabel 'particle values'" << std::endl;
    file << "set ylabel 'frequency'" << std::endl;
    file << std::endl;
    file << "n     = 40" << std::endl;
    file << "max   = "<< v_max << std::endl;
    file << "min   = 0" << std::endl;
    file << "width = (max-min)/(2*n)" << std::endl;
    file << "# function used to map a value to the intervals" << std::endl;
    file << "hist(x,width) = width * floor(x/width) + width/2.0" << std::endl;
    file << std::endl;
    file << "set xrange [min:max]" << std::endl;
    file << "set yrange [0:]" << std::endl;
    file << "set boxwidth 0.5*width" << std::endl;
    file << "set style fill transparent solid 0.4" << std::endl;
    file << "set grid, xtics, ytics" << std::endl;
    file << "plot '"<< data_file << "' using (hist($1,width)):(1.0)";
    file << " lc rgb '#800080' smooth freq with boxes notitle" << std::endl;
  }
  file << std::endl;
  file << "unset multiplot" << std::endl;
  file.close();

  std::string cmd = "gnuplot " + path;
  int exit_code = std::system(cmd.data());
  assert(exit_code == 0);
  std::cout << "Data distribution plots generated in " << output_eps << std::endl;
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

  // parse params
  file >> json;

  input_hacc = json["input"]["filename"];

  for (auto&& scalar : json["input"]["scalars"]) {
    attributes.push_back(scalar);
  }

  num_bins = json["input"]["num_bins"];
  assert(num_bins > 0);

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

    debug_log << "\nLoading and running " << scalar << std::endl;

    if (ioMgr->loadData(scalar)) {
      // save infos for debug
      debug_log << ioMgr->getDataInfo();
      debug_log << ioMgr->getLog();

      size[scalar] = ioMgr->getNumElements();
      computeFrequencies(scalar, ioMgr->data, ioMgr->data);
      computeShannonEntropy(scalar);

    } else {
      if (my_rank == 0)
        std::cerr << "Error while loading " << scalar << ", exiting now" << std::endl;
      return false;
    }
    MPI_Barrier(comm);
  }

  // output logs
  std::string prefix = json["analysis"]["output"]["logs"];
  std::string output = prefix + "_rank_" + std::to_string(my_rank) +".log";

  std::ofstream logfile(output, std::ios::out);
  logfile << debug_log.str();
  debug_log.str("");
  logfile.close();
  std::cout << "Logs generated in "<< output << std::endl;
  MPI_Barrier(comm);

  // dump data and generate plot
  prefix = json["analysis"]["output"]["gnuplot"];
  exportResults(prefix);

  // everything was fine at this point  
  return true;
}

/* -------------------------------------------------------------------------- */


#endif

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
#include "tools.h"
// data loader
#include "dataLoaderInterface.hpp"
#include "HACCDataLoader.hpp"

/* -------------------------------------------------------------------------- */
class Analyzer {

public: 
  Analyzer() = default; 
  Analyzer(Analyzer const&) = delete; 
  Analyzer(Analyzer&&) noexcept = delete; 
  Analyzer(const char* in_path, int in_rank, int in_nb_ranks, MPI_Comm in_comm);
  ~Analyzer() = default;

  bool extractNonHalos(std::string output);
  bool run();  // ok

private:
  bool computeFrequencies(std::string scalar, void* original, void* approx);
  double computeShannonEntropy(std::string scalar);
  void generateHistogram(std::string path = ".");

  // IO
  std::string json_path;
  std::string input_full;
  std::string input_halo;
  std::string output_log;
  std::string output_gnu;
  std::unique_ptr<DataLoaderInterface> ioMgr;
  std::stringstream debug_log;

  // per halo attribute data
  size_t num_bins = 0;
  std::vector<std::string> attributes;
  std::unordered_map<std::string, std::vector<double>> frequency;
  std::unordered_map<std::string, size_t> size;

  // MPI
  int my_rank  = 0;
  int nb_ranks = 0;
  MPI_Comm comm = MPI_COMM_NULL;
};

/* -------------------------------------------------------------------------- */
inline Analyzer::Analyzer(const char* in_path, int in_rank,
                                int in_nb_ranks, MPI_Comm in_comm)
  : json_path(in_path), my_rank(in_rank),
    nb_ranks(in_nb_ranks), comm(in_comm)
{
  assert(nb_ranks > 0);
  // parse parameters
  nlohmann::json json;

  std::ifstream file(json_path);
  assert(file.is_open());
  assert(file.good());

  // parse params
  file >> json;

  if (json["input"].find("halo") != json["input"].end())
    input_halo = json["input"]["halo"];

  if (json["input"].find("full") != json["input"].end())
    input_full = json["input"]["full"];

  for (auto&& scalar : json["input"]["scalars"])
    attributes.push_back(scalar);

  if (json["input"].find("num_bins") != json["input"].end())
    num_bins = json["input"]["num_bins"];
  assert(num_bins > 0);

  // set the IO manager
  ioMgr = std::make_unique<HACCDataLoader>();

  if (json["input"].find("datainfo") != json["input"].end()) {
    auto const& datainfo = json["input"]["datainfo"];
    for (auto it = datainfo.begin(); it != datainfo.end(); ++it) {
      ioMgr->loaderParams[it.key()] = strConvert::toStr(it.value());
    }
  }

  // set outputs paths
  std::string prefix = json["analysis"]["output"]["logs"];
  output_log = prefix + "_rank_" + std::to_string(my_rank) +".log";
  output_gnu = json["analysis"]["output"]["gnuplot"];

}

/* -------------------------------------------------------------------------- */
inline bool Analyzer::computeFrequencies(std::string scalar, void* original, void* approx) {

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
inline double Analyzer::computeShannonEntropy(std::string scalar) {

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
inline void Analyzer::generateHistogram(std::string root_path) {

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
  file << num_bins << " bins - file: " << tools::base(input_halo) << "'" << std::endl;

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
inline bool Analyzer::run() {

  debug_log.clear();
  debug_log.str("");
  debug_log << "Found "<< attributes.size() <<" attributes" << std::endl;

  ioMgr->init(input_halo, comm); 
  ioMgr->setSave(false); 
  MPI_Barrier(comm);

  for (auto&& scalar: attributes) {

    debug_log << "\nLoading and running " << scalar << std::endl;
    // load current data
    assert(ioMgr->loadData(scalar));

    // save infos for debug
    debug_log << ioMgr->getDataInfo();
    debug_log << ioMgr->getLog();

    size[scalar] = ioMgr->getNumElements();
    computeFrequencies(scalar, ioMgr->data, ioMgr->data);
    computeShannonEntropy(scalar);

    MPI_Barrier(comm);
  }

  // output logs
  std::ofstream logfile(output_log, std::ios::out);
  logfile << debug_log.str();
  logfile.close();
  std::cout << "Logs generated in "<< output_log << std::endl;
  MPI_Barrier(comm);

  // dump data and generate plot
  generateHistogram(output_gnu);

  // everything was fine at this point  
  return true;
}

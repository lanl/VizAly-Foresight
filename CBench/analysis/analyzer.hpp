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
#include <set>
#include "tools.h"
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

  std::vector<float> extractNonHalos(std::string scalar, bool debug_dump = false);
  bool run();  // ok

private:
  bool computeFrequencies(std::string scalar, float* original, float* approx);
  double computeShannonEntropy(std::string scalar);
  void generateHistogram();
  void saveDumpInfos(std::string input, int id);
  void dumpLogs();

  // IO
  std::string json_path;
  std::string input_full;
  std::string input_halo;
  std::string output_log;
  std::string output_gnu;
  std::stringstream debug_log;
  std::unique_ptr<HACCDataLoader> ioMgr;

  // per halo attribute data
  bool extract_non_halos = false;
  size_t count_halos = 0;
  size_t count_non_halos = 0;

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

  // check if non halo particles should be extracted
  if (json["analysis"].find("extract_non_halos") != json["analysis"].end()) {
    assert(input_full != "");
    extract_non_halos = (json["analysis"]["extract_non_halos"] == "yes");
  }

  // set outputs paths
  std::string prefix = json["analysis"]["output"]["logs"];
  output_log = prefix + "_rank_" + std::to_string(my_rank) +".log";
  output_gnu = json["analysis"]["output"]["gnuplot"];

}

/* -------------------------------------------------------------------------- */
inline bool Analyzer::computeFrequencies(std::string scalar, float* original, float* approx) {

  double const& n = size[scalar];
  assert(n > 0.);

  // step 1. determine lower and upper bounds on data
  double global_max = 0;
  double global_min = 0;
  double local_max = std::numeric_limits<float>::min();
  double local_min = std::numeric_limits<float>::max();

  for (size_t i = 0; i < n; ++i) {
    if (original[i] > local_max) { local_max = original[i]; }
    if (original[i] < local_min) { local_min = original[i]; }
  }

  MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, comm);

  debug_log << "= local_extents: [" << local_min << ", " << local_max << "]"<< std::endl;
  debug_log << "= global_extents: [" << global_min << ", " << global_max << "]"<< std::endl;
  MPI_Barrier(comm);

  // Compute histogram of values
  if (global_max > 0) {

    debug_log << "num_bins: " << num_bins << std::endl;

    size_t local_histogram[num_bins];
    size_t global_histogram[num_bins];

    std::fill(local_histogram, local_histogram + num_bins, 0);
    std::fill(global_histogram, global_histogram + num_bins, 0);

    double range = global_max - global_min;
    double capacity = range / num_bins;

    for (size_t i = 0; i < n; ++i) {

      double value = approx[i];
      double relative_value = (value - global_min) / range;
      int index = (range * relative_value) / capacity;

      if (index >= num_bins) { index--; }

      local_histogram[index]++;
    }

    MPI_Allreduce(local_histogram, global_histogram, num_bins, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

    // fill frequency eventually
    frequency[scalar].clear();
    frequency[scalar].resize(num_bins);

    for (size_t i = 0; i < num_bins; ++i)
      frequency[scalar][i] = static_cast<double>(global_histogram[i]) / n;

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
inline void Analyzer::generateHistogram() {

  if (my_rank != 0)
    return;

  auto const& root_path = output_gnu;
  auto const num_bins_str = std::to_string(num_bins);

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

  std::string title;
  if (not extract_non_halos) {
    title = "halo data - "+ num_bins_str +" bins - file: "
            + tools::base(input_halo)
            + "("+ std::to_string(count_halos) +" particles)";
  } else {
    title = "non halo data - "+ num_bins_str +" bins - file: "
            + tools::base(input_full)
            + "("+ std::to_string(count_non_halos) +" particles)";
  }

  file << "#" << std::endl;
  file << "# HACC data distribution gnuplot script" << std::endl;
  file << "#" << std::endl;

  file << "reset" << std::endl;
  file << "set terminal postscript eps enhanced color 14 size 21cm,14cm" << std::endl;
  file << "set output '"<< output_eps <<"'" << std::endl;
  file << std::endl;
  file << "set multiplot layout 2,3 title '"<< title <<"'" << std::endl;

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

  if (not extract_non_halos) {

    ioMgr->init(input_halo, comm);
    ioMgr->setSave(false);
    MPI_Barrier(comm);

    for (auto&& scalar: attributes) {
      debug_log << "\nLoading and running " << scalar << std::endl;
      // load current data
      if (ioMgr->loadData(scalar)) {
        // save infos for debug
        debug_log << ioMgr->getDataInfo();
        debug_log << ioMgr->getLog();

        size[scalar] = ioMgr->getNumElements();
        float* data = static_cast<float*>(ioMgr->data);
        computeFrequencies(scalar, data, data);
        computeShannonEntropy(scalar);
        count_halos = size[scalar];
        MPI_Barrier(comm);
      }
    }
  } else {
    count_non_halos = 0;

    for (auto&& scalar: attributes) {
      // first extract non halos data
      auto non_halos = extractNonHalos(scalar);
      size[scalar] = ioMgr->getNumElements();
      computeFrequencies(scalar, non_halos.data(), non_halos.data());
      computeShannonEntropy(scalar);
      non_halos.clear();
      MPI_Barrier(comm);
    }
    count_non_halos /= attributes.size();
  }

  // dump data and generate plot
  generateHistogram();

  return true;
}

/* -------------------------------------------------------------------------- */
inline void Analyzer::dumpLogs() {

  std::ofstream logfile(output_log, std::ios::out);
  logfile << debug_log.str();
  logfile.close();
  std::cout << "Logs generated in "<< output_log << std::endl;

  debug_log.clear();
  debug_log.str("");
  MPI_Barrier(comm);
}
/* -------------------------------------------------------------------------- */
inline void Analyzer::saveDumpInfos(std::string input, int id) {

  gio::GenericIO gioReader(comm, input);
  std::vector<gio::GenericIO::VariableInfo> scalar_info;

  gioReader.openAndReadHeader(gio::GenericIO::MismatchRedistribute);
  assert(gioReader.readNRanks() >= nb_ranks);

  gioReader.readPhysOrigin(ioMgr->physOrigin);
  gioReader.readPhysScale(ioMgr->physScale);
  gioReader.getVariableInfo(scalar_info);

  ioMgr->inOutData.emplace_back(id,
    scalar_info[id].Name,
    static_cast<int>(scalar_info[id].Size),
    scalar_info[id].IsFloat,
    scalar_info[id].IsSigned,
    scalar_info[id].IsPhysCoordX,
    scalar_info[id].IsPhysCoordY,
    scalar_info[id].IsPhysCoordZ
  );

  ioMgr->inOutData.back().determineDataType();

}
/* -------------------------------------------------------------------------- */
// Warning: very memory-consuming routine (may segfault on small nodes).
// It may need redistribution due to a possible data partitions mismatch.
inline std::vector<float> Analyzer::extractNonHalos(std::string scalar, bool debug_dump) {

  debug_log.clear();
  debug_log.str("");

  std::vector<float> all;
  std::vector<float> halo;
  std::vector<float> non_halo;

  // load halo only particles data
  debug_log << "Handling halo particles data ... " << std::endl;

  ioMgr->init(input_halo, comm);

  if (ioMgr->loadData(scalar)) {
    debug_log << ioMgr->getDataInfo();
    size_t const n = ioMgr->getNumElements();
    float* data = static_cast<float*>(ioMgr->data);

    halo.resize(n);
    std::copy(data, data+n, halo.begin());
    std::sort(halo.begin(), halo.end());
    ioMgr->close();

    if (my_rank == 0)
      std::cout << "scalar "<< scalar << " copied and sorted" << std::endl;
  }

  MPI_Barrier(comm);

  // load all particles data
  debug_log << "Handling full particles data ... " << std::endl;

  ioMgr->init(input_full, comm);
  if (debug_dump) {
    ioMgr->setSave(true);
    saveDumpInfos(input_full, 0);
  }

  if (ioMgr->loadData(scalar)) {
    debug_log << ioMgr->getDataInfo();
    size_t const n = ioMgr->getNumElements();
    float* data = static_cast<float*>(ioMgr->data);

    all.resize(n);
    std::copy(data, data + n, all.begin());
    std::sort(all.begin(), all.end());
    ioMgr->close();

    if (my_rank == 0)
      std::cout << "scalar " << scalar << " copied and sorted" << std::endl;
  }

  MPI_Barrier(comm);

  // then take the difference
  non_halo.resize(all.size());
  int const ratio = 100. * (double) halo.size() / all.size();

  debug_log << "compute set difference for scalar '" << scalar << "': "
            << "total particles: " << all.size() << ", "
            << "halo particles: " << halo.size() << " ["<< ratio <<" %]."
            << std::endl;

  auto first = non_halo.begin();
  auto last = std::set_difference(all.begin(), all.end(),
                                  halo.begin(), halo.end(), first);

  unsigned long global_count = 0;
  unsigned long local_count = last - first;
  non_halo.resize(local_count);

  // retrieve number of non halos particles
  MPI_Allreduce(&local_count, &global_count, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);

  debug_log << "= local non halo particles extracted: "<< local_count << std::endl;
  debug_log << "= total non halo particles extracted: "<< global_count << std::endl;

  ioMgr->numElements = local_count;
  ioMgr->totalNumberOfElements = global_count;
  count_non_halos += global_count;

  if (my_rank == 0)
    std::cout << debug_log.str();

  if (debug_dump) {

    all.clear();
    all.shrink_to_fit();
    halo.clear();
    halo.shrink_to_fit();

    // dump data
    debug_log << "Dumping non halo particles data ... " << std::endl;
    if (my_rank == 0)
      std::cout << "Dumping non halo particles data" << std::endl;

    ioMgr->saveCompData(scalar, static_cast<void*>(non_halo.data()));
    ioMgr->writeData("non_halo");
    debug_log << " done" << std::endl;
  }

  dumpLogs();
  return std::move(non_halo);
}
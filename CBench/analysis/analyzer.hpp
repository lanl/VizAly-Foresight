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
#define STORE_PARTICLE_MASK 0
/* -------------------------------------------------------------------------- */
class Analyzer {

public: 
  Analyzer() = default; 
  Analyzer(const char* in_path, int in_rank, int in_nb_ranks, MPI_Comm in_comm);
  Analyzer(Analyzer const&) = delete;
  Analyzer(Analyzer&&) noexcept = delete;
  ~Analyzer() = default;

  bool run();

private:
  bool computeFrequencies(int i, float* data);
  void computeShannonEntropy(int i);
  void filterParticles();
  void extractNonHalos(int i);
  void dumpNonHalosData();
  void dumpLogs();
  void generateHistogram();

  // IO
  std::string json_path;
  std::string input_full;
  std::string input_halo;
  std::string output_log;
  std::string output_gnu;
  std::string output_non_halos;
  std::stringstream debug_log;
  std::unique_ptr<HACCDataLoader> ioMgr;

  int num_scalars = 0;
  int num_bins = 0;
  bool extract_non_halos = false;

  long local_parts = 0;
  long total_parts = 0;
  long local_halos = 0;
  long total_halos = 0;
  long local_non_halos = 0;
  long total_non_halos = 0;
  long local_unmatched = 0;
  long total_unmatched = 0;

  // per-scalar data
  std::vector<std::string> scalars;
  std::vector<size_t> count;
  std::vector<double> entropy;
  std::vector<std::vector<double>> frequency;
  std::vector<std::vector<float>> non_halos;

  // per-particle data
  std::vector<bool> is_halo;
  std::vector<long> non_halos_id;
  std::vector<short> non_halos_mask;

  // MPI
  int my_rank  = 0;
  int nb_ranks = 0;
  MPI_Comm comm = MPI_COMM_NULL;
};

/* -------------------------------------------------------------------------- */
inline Analyzer::Analyzer(const char* in_path, int in_rank, int in_nb_ranks, MPI_Comm in_comm)
  : json_path(in_path),
    my_rank(in_rank),
    nb_ranks(in_nb_ranks),
    comm(in_comm)
{
  assert(nb_ranks > 0);
  nlohmann::json json;
  std::string buffer;

  std::ifstream file(json_path);
  assert(file.is_open());
  assert(file.good());

  // parse params and do basic checks
  file >> json;

  assert(json["input"].find("halo") != json["input"].end());
  assert(json["input"].find("full") != json["input"].end());
  assert(json["input"].find("scalars") != json["input"].end());
  assert(json["analysis"].find("entropy") != json["analysis"].end());
  assert(json["analysis"].find("non-halos") != json["analysis"].end());

  // set the IO manager
  ioMgr = std::make_unique<HACCDataLoader>();
  input_halo = json["input"]["halo"];
  input_full = json["input"]["full"];

  for (auto&& name : json["input"]["scalars"])
    scalars.push_back(name);

  // init data structures
  num_bins = json["analysis"]["entropy"]["num_bins"];
  num_scalars = scalars.size();
  entropy.resize(num_scalars);
  count.resize(num_scalars);
  frequency.resize(num_scalars);

  // set outputs paths
  std::string prefix = json["analysis"]["entropy"]["logs"];
  output_log = prefix + "_rank_" + std::to_string(my_rank) +".log";
  output_gnu = json["analysis"]["entropy"]["plots"];

  // non-halos particles
  // check if non halo particles should be extracted
  std::string option = json["analysis"]["non-halos"]["extract"];

  if (option == "yes") {
    extract_non_halos = true;
    output_non_halos = json["analysis"]["non-halos"]["output"];
    non_halos.resize(num_scalars);
  }
}

/* -------------------------------------------------------------------------- */
inline bool Analyzer::computeFrequencies(int i, float* data) {

  auto const n = (extract_non_halos ? local_non_halos : local_halos);
  assert(n > 0);

  debug_log << "computing frequencies for scalar '"<< scalars[i] <<"'"<< std::endl;

  // step 1. determine lower and upper bounds on data
  double local_max = std::numeric_limits<float>::min();
  double local_min = std::numeric_limits<float>::max();
  double total_max = 0;
  double total_min = 0;

  for (auto j=0; j < n; ++j) {
    if (data[j] > local_max) { local_max = data[j]; }
    if (data[j] < local_min) { local_min = data[j]; }
  }

  MPI_Allreduce(&local_max, &total_max, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(&local_min, &total_min, 1, MPI_DOUBLE, MPI_MIN, comm);

  debug_log << "= local_extents: [" << local_min << ", " << local_max << "]"<< std::endl;
  debug_log << "= total_extents: [" << total_min << ", " << total_max << "]"<< std::endl;
  MPI_Barrier(comm);

  // Compute histogram of values
  if (total_max > 0) {

    debug_log << "num_bins: " << num_bins << std::endl;

    long local_histogram[num_bins];
    long total_histogram[num_bins];

    std::fill(local_histogram, local_histogram + num_bins, 0);
    std::fill(total_histogram, total_histogram + num_bins, 0);

    double range = total_max - total_min;
    double capacity = range / num_bins;

    for (auto k=0; k < n; ++k) {
      double relative_value = (data[k] - total_min) / range;
      int index = (range * relative_value) / capacity;

      if (index >= num_bins)
        index--;

      local_histogram[index]++;
    }

    MPI_Allreduce(local_histogram, total_histogram, num_bins, MPI_LONG, MPI_SUM, comm);

    // fill frequency eventually
    frequency[i].clear();
    frequency[i].resize(num_bins);

    for (int j=0; j < num_bins; ++j)
      frequency[i][j] = double(total_histogram[j]) / n;

    return true;
  }

  return false;
}

/* -------------------------------------------------------------------------- */
inline void Analyzer::computeShannonEntropy(int i) {

  double const epsilon = std::numeric_limits<double>::epsilon();

  entropy[i] = 0.;
  for (auto&& value : frequency[i]) {
    auto const& p_j = std::max(value, epsilon);
    entropy[i] += (p_j * std::log2(p_j));
  }
  entropy[i] *= -1;

  debug_log << "shannon entropy for '"<< scalars[i] << "' = "
            << entropy[i] << " using " << num_bins << " bins." << std::endl;
}

/* -------------------------------------------------------------------------- */
inline void Analyzer::filterParticles() {

  local_non_halos = 0;
  local_unmatched = 0;
  total_non_halos = 0;
  total_unmatched = 0;

  // map from particle id to dataset index
  std::unordered_map<long, long> index;

  debug_log.clear();
  debug_log.str("");
  debug_log << "step 1/4: set particles index mapping ... " << std::endl;

  ioMgr->init(input_full, comm);

  if (ioMgr->loadData("id")) {
    debug_log << ioMgr->getDataInfo();
    // update particles count
    total_parts = ioMgr->totalNumberOfElements;
    local_parts = ioMgr->getNumElements();
    auto data = static_cast<long*>(ioMgr->data);
    // keep track of particles ID to dataset index such that
    // we can easily identify non halo particles when we load
    // another scalar data from the full particles dataset.
    for (auto i=0; i < local_parts; ++i)
      index[data[i]] = i;

    ioMgr->close();

    // init lookup table
    is_halo.resize(local_parts, false);

    debug_log << "= local: "<< local_parts << std::endl;
    debug_log << "= total: "<< total_parts << std::endl << std::endl;
    MPI_Barrier(comm);
  }

  debug_log << "step 2/4: update lookup table ... " << std::endl;

  ioMgr->init(input_halo, comm);

  if (ioMgr->loadData("id")) {
    debug_log << ioMgr->getDataInfo();
    // update particles count
    total_halos = ioMgr->totalNumberOfElements;
    local_halos = ioMgr->getNumElements();
    auto data = static_cast<long*>(ioMgr->data);

    // update halo particles lookup table
    // while counting unmatched particles due to
    // partition mismatch.
    for (auto i=0; i < local_halos; ++i) {
      if (not index.count(data[i])) {
        ++local_unmatched;
      } else {
        auto const& k = index[data[i]];
        is_halo[k] = true;
      }
    }

    ioMgr->close();
    debug_log << "= local: "<< local_halos << std::endl;
    debug_log << "= total: "<< total_halos << std::endl << std::endl;
    MPI_Barrier(comm);
  }

  // early release to reduce memory pressure
  index.clear();

  debug_log << "step 3/4: store non-halos id and mask ... " << std::endl;

  ioMgr->init(input_full, comm);
  ioMgr->setSave(true);
  ioMgr->saveInputFileParameters();

  non_halos_id.reserve(local_non_halos);
  non_halos_mask.reserve(local_non_halos);

  if (ioMgr->loadData("id")) {
    auto* data = static_cast<long*>(ioMgr->data);
    for (auto i=0; i < local_parts; ++i) {
      if (not is_halo[i]) {
        non_halos_id.push_back(data[i]);
        ++local_non_halos;
      }
    }
    ioMgr->close();
  }

#if STORE_PARTICLE_MASK
  MPI_Barrier(comm);

  if (ioMgr->loadData("mask")) {
    auto* data = static_cast<short*>(ioMgr->data);
    for (auto i=0; i < local_parts; ++i) {
      if (not is_halo[i])
        non_halos_mask.push_back(data[i]);
    }
    ioMgr->close();
  }
#endif

  debug_log << "= copied: "<< non_halos_id.size() << std::endl << std::endl;
  MPI_Barrier(comm);

  debug_log << "step 4/4: assess particles count ... " << std::endl;

  // retrieve number of non halos particles
  MPI_Allreduce(&local_non_halos, &total_non_halos, 1, MPI_LONG, MPI_SUM, comm);
  MPI_Allreduce(&local_unmatched, &total_unmatched, 1, MPI_LONG, MPI_SUM, comm);

  double const ratio[] = {
    100. * double(total_halos) / total_parts,
    100. * double(total_non_halos) / total_parts,
    100. * double(total_unmatched) / total_halos
  };

  debug_log << "= total particles: "<< total_parts     << std::endl;
  debug_log << "= total halos:     "<< total_halos     << " ["<< ratio[0] <<" %]." << std::endl;
  debug_log << "= total non-halos: "<< total_non_halos << " ["<< ratio[1] <<" %]." << std::endl;
  debug_log << "= total unmatched: "<< total_unmatched << " ["<< ratio[2] <<" %]." << std::endl;
  debug_log << std::endl;

  if (my_rank == 0)
    std::cout << debug_log.str();

  MPI_Barrier(comm);
  dumpLogs();
}

/* -------------------------------------------------------------------------- */
inline void Analyzer::extractNonHalos(int i) {

  assert(i < num_scalars);

  debug_log << "Cache non-halos particles data '"<< scalars[i] <<"'"<< std::endl;

  ioMgr->filename = input_full;

  if (ioMgr->loadData(scalars[i])) {
    size_t const n = ioMgr->getNumElements();
    assert(n == local_parts);
    non_halos[i].reserve(n);
    auto data = static_cast<float*>(ioMgr->data);

    for (auto k=0; k < local_parts; ++k) {
      if (not is_halo[k]) {
        non_halos[i].push_back(data[k]);
      }
    }

    ioMgr->close();

    debug_log << "= local: "<< n << std::endl;
    debug_log << "= total: "<< total_parts << std::endl << std::endl;
    MPI_Barrier(comm);
  }
}

/* -------------------------------------------------------------------------- */
void Analyzer::dumpNonHalosData() {

  debug_log << "Dumping non halos data into '"<< output_non_halos <<"' ... ";

  int periods[3] = {0,0,0};
  auto dim_size = ioMgr->mpiCartPartitions;
  MPI_Cart_create(comm, 3, dim_size, periods, 0, &comm);

  // init writer and open file
  gio::GenericIO gioWriter(comm, output_non_halos);
  gioWriter.setNumElems(local_non_halos);

  // init physical params
  for (int d=0; d < 3; ++d) {
    gioWriter.setPhysOrigin(ioMgr->physOrigin[d], d);
    gioWriter.setPhysScale(ioMgr->physScale[d], d);
  }

  MPI_Barrier(comm);

  unsigned default_flag = gio::GenericIO::VarHasExtraSpace;

  // populate params now
  for (int i=0; i < num_scalars; ++i) {
    auto scalar = scalars[i].c_str();
    unsigned flag = default_flag;
    switch (i) {
      case 0: flag |= gio::GenericIO::VarIsPhysCoordX; break;
      case 1: flag |= gio::GenericIO::VarIsPhysCoordY; break;
      case 2: flag |= gio::GenericIO::VarIsPhysCoordZ; break;
      default: break;
    }
    gioWriter.addVariable(scalar, non_halos[i].data(), flag);
  }

  gioWriter.addVariable("id", non_halos_id.data(), default_flag);
#if STORE_PARTICLE_MASK
  gioWriter.addVariable("mask", non_halos_mask.data(), default_flag);
#endif
  gioWriter.write();

  // release memory eventually
  non_halos.clear();
  non_halos_id.clear();
  non_halos_mask.clear();

  debug_log << " done." << std::endl;
  MPI_Barrier(comm);
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
inline void Analyzer::generateHistogram() {

  auto const& root_path = output_gnu;
  auto const num_bins_str = std::to_string(num_bins);

  // dump data first
  for (int i = 0; i < num_scalars; ++i) {
    auto const& scalar = scalars[i];
    std::string path = root_path + "_" + scalar + "_" + num_bins_str +".dat";
    std::ofstream file(path, std::ios::out|std::ios::trunc);
    assert(file.is_open());
    assert(file.good());

    file << "# scalar: " << scalar << std::endl;
    file << "# num_bins: " << num_bins_str << std::endl;
    for (auto&& value : frequency[i]) {
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
            + "("+ std::to_string(local_halos) +" particles)";
  } else {
    title = "non halo data - "+ num_bins_str +" bins - file: "
            + tools::base(input_full)
            + "("+ std::to_string(local_non_halos) +" particles)";
  }

  file << "#" << std::endl;
  file << "# HACC data distribution gnuplot script" << std::endl;
  file << "#" << std::endl;

  file << "reset" << std::endl;
  file << "set terminal postscript eps enhanced color 14 size 21cm,14cm" << std::endl;
  file << "set output '"<< output_eps <<"'" << std::endl;
  file << std::endl;
  file << "set multiplot layout 2,3 title '"<< title <<"'" << std::endl;

  for (int i = 0; i < num_scalars; ++i) {
    // get bounds
    double v_max = 0.;
    for (auto&& value : frequency[i]) {
      if (value > v_max) { v_max = value; }
    }

    std::string data_file = root_path + "_" + scalars[i] + "_" + num_bins_str +".dat";

    file << std::endl;
    file << "# ----------------------------------------------------" << std::endl;
    file << "set title '"<< scalars[i] << ", entropy=" << entropy[i] <<"'"<< std::endl;
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
  debug_log << "Found "<< num_scalars <<" attributes" << std::endl;

  if (not extract_non_halos) {

    ioMgr->init(input_halo, comm);
    ioMgr->setSave(false);

    for (int i = 0; i < num_scalars; ++i) {
      auto const& scalar = scalars[i];
      debug_log << "\nLoading and running " << scalar << std::endl;
      // load current data
      if (ioMgr->loadData(scalar)) {
        // save infos for debug
        debug_log << ioMgr->getDataInfo();
        debug_log << ioMgr->getLog();
        // set halo particles count
        if (i == 0) {
          local_halos = ioMgr->getNumElements();
          MPI_Allreduce(&local_halos, &total_halos, 1, MPI_LONG, MPI_SUM, comm);
        }

        auto data = static_cast<float*>(ioMgr->data);
        computeFrequencies(i, data);
        computeShannonEntropy(i);
        MPI_Barrier(comm);
      }
    }
  } else {
    // set non halos lookup table and store metadata
    filterParticles();
    for (int i = 0; i < num_scalars; ++i) {
      // extract and store non-halo scalar data
      extractNonHalos(i);
      computeFrequencies(i, non_halos[i].data());
      computeShannonEntropy(i);
      MPI_Barrier(comm);
    }

    dumpNonHalosData();
  }

  // dump data and generate plot
  if (my_rank == 0)
    generateHistogram();

  return true;
}
/* -------------------------------------------------------------------------- */
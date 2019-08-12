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
  Analyzer(const char* in_path, int in_rank, int in_nb_ranks, MPI_Comm in_comm);
  Analyzer(Analyzer const&) = delete;
  Analyzer(Analyzer&&) noexcept = delete;
  ~Analyzer();

  bool run();

private:
  bool computeFrequencies(int id, float* original, float* approx);
  void computeShannonEntropy(int i);
  void updateParticlesCount();
  // non-halo particles
  void extractNonHalos(int i);
  void dumpNonHalosData();
  void generateHistogram();
  void dumpLogs();

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
  size_t num_bins = 0;
  double phys_origin[3];
  double phys_scale[3];

  size_t count_halos = 0;
  size_t count_non_halos = 0;
  bool extract_non_halos = false;

  // per scalar data
  std::vector<std::string> attributes;
  std::vector<size_t> count;
  std::vector<std::vector<double>> frequency;
  std::vector<double> entropy;
  std::vector<float*> non_halos;

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

  for (auto&& scalar : json["input"]["scalars"])
    attributes.push_back(scalar);

  // init data structures
  num_bins = json["analysis"]["entropy"]["num_bins"];
  num_scalars = attributes.size();
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
    for (auto&& data : non_halos)
      data = nullptr;
  }
}

/* -------------------------------------------------------------------------- */
inline Analyzer::~Analyzer() {
  for (auto&& data : non_halos)
    delete [] data;
}

/* -------------------------------------------------------------------------- */
inline void Analyzer::updateParticlesCount() {

  // non-halo particles count may slightly differ due to partition mismatch
  // need to take the minimum to avoid issues when dumping them.
  if (extract_non_halos) {
    std::accumulate(count.begin(), count.end(), count_non_halos,
                    [](auto& a, auto& b) { return std::min(a,b); });
  } else {
    count_halos = count[0];
  }
}

/* -------------------------------------------------------------------------- */
inline bool Analyzer::computeFrequencies(int id, float* original, float* approx) {

  double const& n = count[id];
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
    frequency[id].clear();
    frequency[id].resize(num_bins);

    for (size_t i = 0; i < num_bins; ++i)
      frequency[id][i] = static_cast<double>(global_histogram[i]) / n;

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

  debug_log << "shannon entropy for '"<< attributes[i] << "' = "
            << entropy[i] << " using " << num_bins << " bins." << std::endl;
}

/* -------------------------------------------------------------------------- */
inline void Analyzer::generateHistogram() {

  if (my_rank != 0)
    return;

  updateParticlesCount();

  auto const& root_path = output_gnu;
  auto const num_bins_str = std::to_string(num_bins);

  // dump data first
  for (int i = 0; i < num_scalars; ++i) {
    auto const& scalar = attributes[i];
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

  for (int i = 0; i < num_scalars; ++i) {
    // get bounds
    double v_max = 0.;
    for (auto&& value : frequency[i]) {
      if (value > v_max) { v_max = value; }
    }

    std::string data_file = root_path + "_" + attributes[i] + "_" + num_bins_str +".dat";

    file << std::endl;
    file << "# ----------------------------------------------------" << std::endl;
    file << "set title '"<< attributes[i] << ", entropy=" << entropy[i] <<"'"<< std::endl;
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
      auto const& scalar = attributes[i];
      debug_log << "\nLoading and running " << scalar << std::endl;
      // load current data
      if (ioMgr->loadData(scalar)) {
        // save infos for debug
        debug_log << ioMgr->getDataInfo();
        debug_log << ioMgr->getLog();
        count[i] = ioMgr->getNumElements();

        float* data = static_cast<float*>(ioMgr->data);
        computeFrequencies(i, data, data);
        computeShannonEntropy(i);
        MPI_Barrier(comm);
      }
    }
  } else {
    for (int i = 0; i < num_scalars; ++i) {
      // first extract non halos data
      extractNonHalos(i);
      count[i] = ioMgr->getNumElements();
      computeFrequencies(i, non_halos[i], non_halos[i]);
      computeShannonEntropy(i);
      MPI_Barrier(comm);
    }

    dumpNonHalosData();
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
// warning: very memory-consuming routine (may segfault on small nodes).
inline void Analyzer::extractNonHalos(int i) {

  assert(i < num_scalars);

  std::vector<float> all;
  std::vector<float> halo;
  std::vector<float> unmatched;

  size_t total_halo_count = 0;
  size_t total_full_count = 0;

  // load halo only particles data
  debug_log.clear();
  debug_log.str("");
  debug_log << "Handling halo particles data ... " << std::endl;

  ioMgr->init(input_halo, comm);
  auto scalar = attributes[i];

  if (ioMgr->loadData(scalar)) {
    debug_log << ioMgr->getDataInfo();
    total_halo_count = ioMgr->totalNumberOfElements;
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
  ioMgr->saveInputFileParameters();
  ioMgr->setSave(true);

  if (ioMgr->loadData(scalar)) {

    if (my_rank == 0) {
      std::cout << ioMgr->getLog() << std::endl;
      std::cout << "===================================" << std::endl;
    }

    debug_log << ioMgr->getDataInfo();
    total_full_count = ioMgr->totalNumberOfElements;
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
  //non_halo.resize(all.size());
  non_halos[i] = new float[all.size()];
  double const ratio = 100. * double(total_halo_count) / total_full_count;

  debug_log << "compute set difference for scalar '" << scalar << "': "
            << "all: " << total_full_count << ", "
            << "halo: "<< total_halo_count << " ["<< ratio <<" %]." << std::endl;

  auto first = non_halos[i];//non_halo.begin();
  auto last = std::set_difference(all.begin(), all.end(),
                                  halo.begin(), halo.end(), first);

  unsigned long total_non_halo_count = 0;
  unsigned long local_non_halo_count = last - first;
  //non_halo.resize(local_non_halo_count);

  // retrieve number of non halos particles
  MPI_Allreduce(&local_non_halo_count, &total_non_halo_count, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);

  // now compute the error
  debug_log << "compute unmatched halo particles "<< std::endl;

  unmatched.resize(halo.size());
  auto first_pos = unmatched.begin();
  auto last_pos  = std::set_difference(halo.begin(), halo.end(),
                                       all.begin(), all.end(), first_pos);

  unsigned long total_unmatched_halo_count = 0;
  unsigned long local_unmatched_halo_count = last_pos - first_pos;

  // retrieve number of non halos particles
  MPI_Allreduce(&local_unmatched_halo_count, &total_unmatched_halo_count, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);
  // deduce error
  double const total_error_ratio = 100. * double(total_unmatched_halo_count) / total_halo_count;
  // ---------------

  ioMgr->numElements = local_non_halo_count;
  ioMgr->totalNumberOfElements = total_non_halo_count;

  debug_log << "= local unmatched halo particles: "<< local_unmatched_halo_count << std::endl;
  debug_log << "= local non halo particles extracted: "<< local_non_halo_count << std::endl;
  debug_log << std::endl;
  debug_log << "= total unmatched halo particles: "<< total_unmatched_halo_count << std::endl;
  debug_log << "= total non halo particles extracted: "<< total_non_halo_count << std::endl;
  debug_log << "= total extraction error ratio: "<< total_error_ratio << " %"<< std::endl;

  if (my_rank == 0)
    std::cout << debug_log.str();

  dumpLogs();
  //non_halos[i] = std::move(non_halo.d);
}

/* --------------------------------------------------------------------------
void Analyzer::cacheNonHalosData(int i) {

  assert(i < num_scalars);
  assert(data != nullptr);
  assert(count[i]);

  debug_log << "Caching non halos data for scalar "<< attributes[i] <<" ... ";
  // nb: local to the current rank
  non_halos[i] = new float[count[i]];
  std::copy(data, data + count[i], non_halos[i]);

  debug_log << count[i] << " particle data stored." << std::endl;
}*/
/* -------------------------------------------------------------------------- */
void Analyzer::dumpNonHalosData() {

  debug_log << "Dumping non halos data into '"<< output_non_halos <<"' ... ";

  int periods[3] = {0,0,0};
  int const* dim_size = ioMgr->mpiCartPartitions;
  MPI_Cart_create(comm, 3, dim_size, periods, 0, &comm);

  // init writer and open file
  gio::GenericIO gioWriter(comm, output_non_halos);
  gioWriter.setNumElems(count_non_halos);

  // init physical params
  for (int d=0; d < 3; ++d) {
    gioWriter.setPhysOrigin(ioMgr->physOrigin[d], d);
    gioWriter.setPhysScale(ioMgr->physScale[d], d);
  }

  MPI_Barrier(comm);

  auto flag = [](int i) {
    unsigned f = gio::GenericIO::VarHasExtraSpace;
         if (i == 0) { f |= gio::GenericIO::VarIsPhysCoordX; }
    else if (i == 1) { f |= gio::GenericIO::VarIsPhysCoordY; }
    else if (i == 2) { f |= gio::GenericIO::VarIsPhysCoordZ; }
    return f;
  };

  // populate params now
  for (int i=0; i < num_scalars; ++i) {
    auto scalar = attributes[i].data();
    gioWriter.addVariable(scalar, non_halos[i], flag(i));
  }

  gioWriter.write();

  debug_log << " done." << std::endl;
}
/* -------------------------------------------------------------------------- */
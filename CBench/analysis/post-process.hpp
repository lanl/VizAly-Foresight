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
#include "HACCDataLoader.hpp"
/* -------------------------------------------------------------------------- */
class Merger {

public:
  Merger() = default;
  Merger(const char* in_path, int in_rank, int in_nb_ranks, MPI_Comm in_comm);
  Merger(Merger const&) = delete;
  Merger(Merger&&) noexcept = delete;
  ~Merger() = default;

  void run();

private:

  size_t cache(std::vector<float>* dataset, std::vector<long>& index);
  void dump();
  void dumpLogs();

  std::string json_path;
  std::string halo_file;
  std::string non_halo_file;
  std::string output_combined;
  std::string output_log;
  std::stringstream debug_log;
  std::unique_ptr<HACCDataLoader> ioMgr;

  int num_scalars = 0;
  long local_parts = 0;
  long total_parts = 0;
  long local_halos = 0;
  long total_halos = 0;
  long local_non_halos = 0;
  long total_non_halos = 0;

  std::vector<std::string> scalars;
  std::vector<std::vector<float>> halos;
  std::vector<std::vector<float>> non_halos;

  std::vector<long> halos_id;
  std::vector<long> non_halos_id;

  // MPI
  int my_rank  = 0;
  int nb_ranks = 0;
  MPI_Comm comm = MPI_COMM_NULL;
};

/* -------------------------------------------------------------------------- */
inline Merger::Merger(const char* in_path, int in_rank, int in_nb_ranks, MPI_Comm in_comm)
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

  assert(json["input"].find("scalars") != json["input"].end());
  assert(json["post-process"].find("halos") != json["post-process"].end());
  assert(json["post-process"].find("non-halos") != json["post-process"].end());
  assert(json["post-process"].find("output") != json["post-process"].end());

  halo_file = json["post-process"]["halos"];
  non_halo_file = json["post-process"]["non-halos"];
  output_combined = json["post-process"]["output"];
  output_log = json["post-process"]["logs"];

  for (auto&& name : json["input"]["scalars"])
    scalars.push_back(name);

  num_scalars = scalars.size();
  halos.resize(num_scalars);
  non_halos.resize(num_scalars);

  // set the IO manager
  ioMgr = std::make_unique<HACCDataLoader>();
}

/* -------------------------------------------------------------------------- */
// generic method for caching a given dataset
inline size_t Merger::cache(std::vector<float>* dataset, std::vector<long>& index) {

  debug_log << "Caching dataset ... ";
  assert(dataset != nullptr);

  for (int i=0; i < num_scalars; ++i) {
    if (ioMgr->loadData(scalars[i])) {
      auto const n = ioMgr->getNumElements();
      auto const data = static_cast<float*>(ioMgr->data);
      dataset[i].resize(n);
      std::copy(data, data + n, dataset[i].data());
    }
    MPI_Barrier(comm);
  }

  if (ioMgr->loadData("id")) {
    auto const n = ioMgr->getNumElements();
    auto const data = static_cast<long*>(ioMgr->data);
    index.resize(n);
    std::copy(data, data + n, index.data());
  }

  debug_log << " done." << std::endl;
  MPI_Barrier(comm);

  return index.size();
};

/* -------------------------------------------------------------------------- */
inline void Merger::dump() {

  debug_log << "Dumping dataset ... ";

  assert(local_parts > 0);

  int periods[3] = {0,0,0};
  auto dim_size = ioMgr->mpiCartPartitions;
  MPI_Cart_create(comm, 3, dim_size, periods, 0, &comm);

  // init writer and open file
  gio::GenericIO gioWriter(comm, output_combined);
  gioWriter.setNumElems(local_parts);

  // init physical params
  for (int d=0; d < 3; ++d) {
    gioWriter.setPhysOrigin(ioMgr->physOrigin[d], d);
    gioWriter.setPhysScale(ioMgr->physScale[d], d);
  }

  MPI_Barrier(comm);

  unsigned default_flag = gio::GenericIO::VarHasExtraSpace;

  // populate params now
  for (int i=0; i < num_scalars; ++i) {
    unsigned flag = default_flag;
    switch (i) {
      case 0: flag |= gio::GenericIO::VarIsPhysCoordX; break;
      case 1: flag |= gio::GenericIO::VarIsPhysCoordY; break;
      case 2: flag |= gio::GenericIO::VarIsPhysCoordZ; break;
      default: break;
    }
    gioWriter.addVariable(scalars[i].data(), halos[i].data(), flag);
    gioWriter.addVariable(scalars[i].data(), non_halos[i].data(), flag);
  }

  gioWriter.addVariable("id", halos_id.data(), default_flag);
  gioWriter.addVariable("id", non_halos_id.data(), default_flag);
  gioWriter.write();

  // release memory eventually
  halos.clear();
  halos_id.clear();
  non_halos.clear();
  non_halos_id.clear();

  debug_log << " done." << std::endl;
  MPI_Barrier(comm);
}


/* -------------------------------------------------------------------------- */
inline void Merger::run() {

  debug_log << "Reconstruct decompressed dataset:" << std::endl;

  ioMgr->init(halo_file, comm);
  local_halos = cache(halos.data(), halos_id);

  ioMgr->init(non_halo_file, comm);
  ioMgr->setSave(true);
  ioMgr->saveInputFileParameters();

  local_non_halos = cache(non_halos.data(), non_halos_id);

  // get total number of particles
  total_halos = 0;
  total_non_halos = 0;
  MPI_Allreduce(&local_halos, &total_halos, 1, MPI_LONG, MPI_SUM, comm);
  MPI_Allreduce(&local_non_halos, &total_non_halos, 1, MPI_LONG, MPI_SUM, comm);

  local_parts = local_halos + local_non_halos;
  total_parts = total_halos + total_non_halos;
  debug_log << "\tlocal: "<< local_parts << " particles"<< std::endl;
  debug_log << "\ttotal: "<< total_parts << " particles"<< std::endl;

  MPI_Barrier(comm);

  // now dump everything
  dump();
  dumpLogs();
}

/* -------------------------------------------------------------------------- */
inline void Merger::dumpLogs() {

  std::ofstream logfile(output_log, std::ios::out);
  logfile << debug_log.str();
  logfile.close();
  std::cout << "Logs generated in "<< output_log << std::endl;

  debug_log.clear();
  debug_log.str("");
  MPI_Barrier(comm);
}
/* -------------------------------------------------------------------------- */


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

// data loader
#include "dataLoaderInterface.hpp"
#include "HACCDataLoader.hpp"
/* -------------------------------------------------------------------------- */
class HaloData {

  // mimic 'uncompressed data' from NYX
}

/* -------------------------------------------------------------------------- */
class HaloEntropy {

  // Workflow
  // - load JSON param file
  // - load HACC data and store attributes
  // - test: show attributes
  // - compute Shannon entropy distribution
  // - generate gnuplot histogram script 
  // (!) make it self-contained

private:
  // hacc data loader
  DataLoaderInterface *ioMgr;
  // a map to store shannon entropy for each hacc attribute
  // fixme: not sure about those attributes, to be updated
  std::unordered_map<std::string, double> entropy;
  
  // MPI
  int my_rank;
  int nb_ranks;
  int thread_support;

public:  
  HaloEntropy();
  ~HaloEntropy();
 
  bool init(int argc, char* argv[]); 
  bool loadData(std::string jsonParamFile);
  bool showAttributeData() const;
  bool computeShannonEntropyDistribution(std::string attribName);
  bool generateGnuplot(std::string scriptFile) const;

};
/* -------------------------------------------------------------------------- */
inline HaloEntropy::HaloEntropy() {
  my_rank = 0;
  nb_ranks = 0;
  ioMgr = new HACCDataLoader();
  entropy.clear();
}
/* -------------------------------------------------------------------------- */
inline HaloEntropy::~HaloEntropy() {
  if (ioMgr != nullptr)
    delete ioMgr;
}
/* -------------------------------------------------------------------------- */
inline HaloEntropy::init(int argc, char* argv[]) {

  // initialize MPI environment
  MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &thread_support);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
}
/* -------------------------------------------------------------------------- */
inline bool HaloEntropy::loadData(std::string path) {

  nlohmann::json json;
  std::ifstream file(path);
  std::stringstream debuglog;
  
  std::string filetype;
  std::string input_hacc;
  std::string output_log;
  std::vector<std::string> attributes;


  try {
    // pass file to json parser
    file >> json;
  } catch(nlohmann::json::parse_error& e) {
    if (my_rank == 0) 
      std::cerr << "Invalid JSON file " << path << "\n" << e.what() << std::endl;
    return false;
  }

  // retrieve input file
  filetype = json["input"]["filetype"];
  input_hacc = json["input"]["filename"];
  if (filetype != "HACC") {
    if (my_rank == 0) 
      std::cerr << "Only HACC data is supported" << std::endl;
    return false;
  }

  // setup logs
  log_file = json["cbench"]["output"]["log-file"]; 
  output_log = "logs/"+ log_file +"_"+ std::to_string(my_rank);
  writeLog(output_log, debuglog.str());

  for (auto&& scalar : json["input"]["scalars"]) {
    attributes.emplace_back(scalar);
  }

  // set the IO manager
  assert(ioMgr != nullptr);

  if (json["input"].find("datainfo") != json["input"].end()) {
    auto const& datainfo = json["input"]["datainfo"];
    for (auto it = datainfo.begin(); it != datainfo.end(); ++it) {
      ioMgr->loaderParams[it.key()] = strConvert::toStr(it.value());
    }    
  }

  ioMgr->init(input_hacc, MPI_COMM_WORLD); 
  ioMgr->setSave(false); 
  MPI_Barrier(MPI_COMM_WORLD);

  // loop over all scalars and load each of them
  for (auto&& scalar: scalars) {
    debuglog << "\nLoading " << scalar << std::endl;
    if (ioMgr->loadData(scalar)) {
      // save infos for debug
      debuglog << ioMgr->getDataInfo();
      debuglog << ioMgr->getLog();
    } else {
      if (my_rank == 0)
        std::cerr << "Error while loading " << scalar << ", exiting now" << std::endl;
      return false;
    }
    appendLog(output_log, debuglog.str());
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // print some metadata on loaded HACC file


  // everything was fine at this point  
  return true;
}

#endif

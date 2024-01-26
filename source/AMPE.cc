// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE.
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
//
//
#include "AMPE.h"
#include "QuatModel.h"

#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/InputManager.h"


AMPE::AMPE(MPI_Comm comm)
{
   tbox::SAMRAI_MPI::init(comm);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
}

AMPE::~AMPE()
{
   delete d_pfm;

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();
}

void AMPE::initialize(const std::string input_filename,
                      const std::string restart_read_dirname,
                      const int restore_num)
{
   std::shared_ptr<tbox::MemoryDatabase> input_db(
       new tbox::MemoryDatabase("input_db"));
   tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

   std::string run_name;
   if (input_db->keyExists("run_name")) {
      run_name = input_db->getString("run_name");
   } else {
      // make from input file name
      size_t a = input_filename.find_last_of("/");
      if (a == input_filename.npos) a = 0;
      if (a > 0) a++;
      size_t b = input_filename.rfind(".");
      run_name = input_filename.substr(a, b - a);
      tbox::pout << "run_name = " << run_name << std::endl;
   }

   bool log_all_nodes = false;
   std::string log_file_name = run_name + ".log";

   if (input_db->isDatabase("Logging")) {
      std::shared_ptr<tbox::Database> log_db = input_db->getDatabase("Logging");
      if (log_db->keyExists("filename")) {
         log_file_name = log_db->getString("filename");
      }

      if (log_db->keyExists("log_all_nodes")) {
         log_all_nodes = log_db->getBool("log_all_nodes");
      }
   }

   if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_file_name);
   } else {
      tbox::PIO::logOnlyNodeZero(log_file_name);
   }

   tbox::plog << "AMPE: git_version " << gitCommitID() << std::endl;

   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());
   tbox::plog << "Run with " << mpi.getSize() << " MPI tasks" << std::endl;

   bool is_from_restart = (restore_num >= 0);
   if (is_from_restart) {
      tbox::plog << "restart_read_dirname = " << restart_read_dirname
                 << std::endl;
      tbox::plog << "restore_num = " << restore_num << std::endl;
   }

   // Create timers
   tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

   d_time_man = tbox::TimerManager::getManager();
   d_time_man->resetAllTimers();

   // Create a PFModel object
   std::string model_type = input_db->getStringWithDefault("model_type",
                                                           "Qua"
                                                           "t");
   if (model_type == "Quat") {
      d_pfm = new QuatModel(4);
   }
#if NDIM == 2
   else if (model_type == "KWC") {
      d_pfm = new QuatModel(1);
   } else if (model_type == "KWCcomplex") {
      d_pfm = new QuatModel(2);
   }
#endif
   else {
      TBOX_ERROR("Invalid model_type" << std::endl);
   }

   d_pfm->Initialize(input_db, run_name, is_from_restart, restart_read_dirname,
                     restore_num);
   /*
    * After creating all objects and initializing their state, we
    * print the input database and variable database contents to
    * the log file.
    */
   tbox::plog << "\nCheck input data and variables before simulation:"
              << std::endl;
   tbox::plog << "Input database..." << std::endl;
   input_db->printClassData(tbox::plog);
   tbox::plog << "\nVariable database..." << std::endl;
   hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

   input_db.reset();
}

void AMPE::run()
{
   d_pfm->Run();

   d_time_man->print(tbox::plog);
}

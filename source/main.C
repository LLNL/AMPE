// Copyright (c) 2018, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
// Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
// LLNL-CODE-747500
// All rights reserved.
// This file is part of AMPE. 
// For details, see https://github.com/LLNL/AMPE
// Please also read AMPE/LICENSE.
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
// LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
// 
#include "SAMRAI/SAMRAI_config.h"

#include <string>
#include <fstream>

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include <boost/make_shared.hpp>
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Database.h"

#include "ran2.h"
#include "PFModel.h"
#include "QuatModel.h"

using namespace SAMRAI;
using namespace std;

int main( int argc, char *argv[] )
{
   // Initialize MPI, SAMRAI, and enable logging.

   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();
   const tbox::SAMRAI_MPI& mpi(tbox::SAMRAI_MPI::getSAMRAIWorld());

   /* This extra code block is used to scope some temporaries that are
    * created, it forces the destruction before the manager is
    * shutdown.
    */
   {

   //-----------------------------------------------------------------------
   /*
    * Process command line arguments and dump to log file.
    * For non-restarted case, command line is:
    *
    *    executable <input file name>
    *
   * For restarted run, command line is:
    *
    *    executable <input file name> <restart directory> \
    *               <restart number>
    */

   string input_filename;
   string restart_read_dirname;
   int restore_num = 0;
   bool is_from_restart = false;

   if ( (argc != 2) && (argc != 4) ) {
      tbox::pout << "USAGE:  " << argv[0] << " <input filename> " <<
         "<restart dir> <restore number> [options]\n" <<
         "  options:\n" <<
         "  none at this time" << endl;
      tbox::SAMRAI_MPI::abort();
      return (-1);
   }
   else {
      input_filename = argv[1];
      if (argc == 4) {
         restart_read_dirname = argv[2];
         restore_num = atoi( argv[3] );

         is_from_restart = true;
      }
   }

   //-----------------------------------------------------------------------
   // Create input database and parse all data in input file.

   boost::shared_ptr<tbox::MemoryDatabase> input_db(new tbox::MemoryDatabase("input_db"));
   tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

   //-----------------------------------------------------------------------
   // Read key input settings
   
   string run_name;
   if ( input_db->keyExists( "run_name" ) ) {
      run_name = input_db->getString( "run_name" );
   }
   else {
      // make from input file name
      run_name = input_filename.substr( 0, input_filename.rfind( "." ) );
   }

   //-----------------------------------------------------------------------
   //
   // Logfile
   //
   bool log_all_nodes = false;
   string log_file_name = run_name + ".log";

   if ( input_db->isDatabase( "Logging" ) ) {
      boost::shared_ptr<tbox::Database> log_db = input_db->getDatabase( "Logging" );

      if ( log_db->keyExists( "filename" ) ) {
         log_file_name = log_db->getString( "filename" );
      }

      if ( log_db->keyExists( "log_all_nodes" ) ) {
         log_all_nodes = log_db->getBool( "log_all_nodes" );
      }
   }

   if ( log_all_nodes ) {
      tbox::PIO::logAllNodes( log_file_name );
   }
   else {
      tbox::PIO::logOnlyNodeZero( log_file_name );
   }

#ifdef GITVERSION
#define xstr(x) #x
#define LOG(x) tbox::plog<<" AMPE: git_version "<<xstr(x)<<endl;
    LOG(GITVERSION);
    tbox::plog<<endl;
#endif

   tbox::plog << "Run with "<<mpi.getSize()<<" MPI tasks"<<endl; 	
   tbox::plog << "input_filename = " << input_filename << endl;

   if ( is_from_restart ) {
      tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
      tbox::plog << "restore_num = " << restore_num << endl;
   }

   //-----------------------------------------------------------------------
   // Create timers

   tbox::TimerManager::createManager(input_db->getDatabase("TimerManager"));

   tbox::TimerManager* time_man = tbox::TimerManager::getManager();
   time_man->resetAllTimers();

   //-----------------------------------------------------------------------
   // Create a PFModel object

   string model_type = input_db->getStringWithDefault( "model_type", "Quat" );

   PFModel* pfm = NULL;
   if ( model_type == "Quat" ) {
      pfm = new QuatModel( 4 );
   }
#if NDIM==2
   else if ( model_type == "KWC" ) {
      //pfm = new KWCModel();
      pfm = new QuatModel( 1 );
   }
   else if ( model_type == "KWCcomplex" ) {
      pfm = new QuatModel( 2 );
   }
#endif
   else {
      TBOX_ERROR( "Invalid model_type" << endl );
   }

   pfm->Initialize(
      input_db, run_name,
      is_from_restart, restart_read_dirname, restore_num
      );

   /*
    * After creating all objects and initializing their state, we
    * print the input database and variable database contents to
    * the log file.
    */
   tbox::plog << "\nCheck input data and variables before simulation:" << endl;
   tbox::plog << "Input database..." << endl;
   input_db->printClassData(tbox::plog);
   tbox::plog << "\nVariable database..." << endl;
   hier::VariableDatabase::getDatabase()->printClassData(tbox::plog);

   pfm->Run();

   input_db.reset();

   delete pfm;

   time_man->print(tbox::plog);
   //time_man->print(tbox::pout);

   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return(0);
}

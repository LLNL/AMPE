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
#include "CALPHADFreeEnergyFunctions.h"
#include "CALPHADMobility.h"
#include "CompositionStrategyMobilities.h"
#include "CALPHADFreeEnergyStrategy.h"
#include "QuatModelParameters.h"
#include "ConstantMolarVolumeStrategy.h"

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include <boost/make_shared.hpp>
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/TimerManager.h"
#include "SAMRAI/tbox/Database.h"

#include <string>
#include <fstream>

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
    */

   std::string input_filename;
   input_filename = argv[1];

   //-----------------------------------------------------------------------
   // Create input database and parse all data in input file.

   boost::shared_ptr<tbox::MemoryDatabase> input_db(new tbox::MemoryDatabase("input_db"));
   tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

   //-----------------------------------------------------------------------
   // Read key input settings
   
   std::string run_name;
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
   std::string log_file_name = run_name + ".log";

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

#ifdef SVNVERSION
#define xstr(x) #x
#define LOG(x) tbox::plog<<" AMPE: svn_version "<<xstr(x)<<endl;
    LOG(SVNVERSION);
    tbox::plog<<endl;
#endif

   tbox::plog << "Run with "<<mpi.getSize()<<" MPI tasks"<<endl; 	
   tbox::plog << "input_filename = " << input_filename << endl;

   boost::shared_ptr<tbox::Database> model_db =
      input_db->getDatabase("ModelParameters");

   double phase_well_scale = model_db->getDouble( "phi_well_scale" );
   double eta_well_scale   = model_db->getDoubleWithDefault( "eta_well_scale", 0. );

   string eta_well_func_type =
      model_db->getStringWithDefault( "eta_well_func_type", "double" );
   string phase_well_func_type =
         model_db->getString( "phi_well_func_type" );
   if ( eta_well_func_type[0] != 's' &&
        eta_well_func_type[0] != 'S' &&
        eta_well_func_type[0] != 'd' &&
        eta_well_func_type[0] != 'D' ) {
      TBOX_ERROR( "Error: invalid value for eta_well_func_type" );
   }

   string phase_interp_func_type = "pbg";
   string eta_interp_func_type   ="pbg";
   
   boost::shared_ptr<tbox::Database> temperature_db = model_db->getDatabase( "Temperature" );
   double temperature = temperature_db->getDouble( "temperature" );

   boost::shared_ptr<tbox::Database> conc_db(model_db->getDatabase( "ConcentrationModel" ));
   string conc_avg_func_type =
      conc_db->getStringWithDefault( "avg_func_type", "a" );

   boost::shared_ptr<tbox::Database> dcalphad_db=conc_db->getDatabase( "Calphad" );
   std::string calphad_filename = dcalphad_db->getString( "filename" );
   boost::shared_ptr<tbox::MemoryDatabase> calphad_db ( new tbox::MemoryDatabase( "calphad_db" ) );
   tbox::InputManager::getManager()->parseInputFile( calphad_filename, calphad_db );
   
   boost::shared_ptr<tbox::Database> newton_db;
   if ( conc_db->isDatabase( "NewtonSolver" ) )
      newton_db = conc_db->getDatabase( "NewtonSolver" );

   bool with_third_phase=false;
   
   CALPHADFreeEnergyFunctions
      cafe(calphad_db, newton_db,
           phase_interp_func_type,
           eta_interp_func_type,
           conc_avg_func_type,
           with_third_phase,
           phase_well_scale,
           eta_well_scale,
           phase_well_func_type,
           eta_well_func_type);
   
   cafe.printEnergyVsComposition(temperature);

   // choose pair of phases: phaseL, phaseA, phaseB
   const PHASE_INDEX pi0=phaseL;
   const PHASE_INDEX pi1=phaseA;
   
   // initial guesses
   double ceq_init0=0.5;
   double ceq_init1=0.5;

   double lceq[2]={ceq_init0,ceq_init1};
   
   // compute equilibrium concentrations
   bool found_ceq =
      cafe.computeCeqT(temperature,pi0,pi1,&lceq[0]);
   if( lceq[0]>1. )found_ceq = false;
   if( lceq[0]<0. )found_ceq = false;
   if( lceq[1]>1. )found_ceq = false;
   if( lceq[1]<0. )found_ceq = false;
   
   if( !found_ceq )
   {
      lceq[0]=ceq_init1;
      lceq[1]=ceq_init0;
      found_ceq =
         cafe.computeCeqT(temperature,pi0,pi1,&lceq[0]);
   }
   
   if( found_ceq ){
      tbox::plog<<"Found equilibrium concentrations: "<<lceq[0]<<" and "<<lceq[1]<<"..."<<endl;
   }else{
      tbox::plog<<"WARNING: Equilibrium concentrations not found... "<<endl;
   }
   
   QuatModelParameters model_parameters;
   model_parameters.readModelParameters(model_db);

   tbox::plog<<"ConstantMolarVolumeStrategy... "<<endl;
   ConstantMolarVolumeStrategy mvstrategy(model_parameters.molar_volume_liquid(),
                                          model_parameters.molar_volume_solid_A(),
                                          model_parameters.molar_volume_solid_B());
   tbox::plog<<"CALPHADFreeEnergyStrategy... "<<endl;
   CALPHADFreeEnergyStrategy free_energy_strategy(
               calphad_db, newton_db,
               model_parameters.phase_interp_func_type(),
               model_parameters.eta_interp_func_type(),
               model_parameters.conc_avg_func_type(),
               &mvstrategy,
               0,
               1,
               2,
               model_parameters.with_third_phase(),
               model_parameters.phase_well_scale(),
               model_parameters.eta_well_scale(),
               model_parameters.phase_well_func_type(),
               model_parameters.eta_well_func_type() );

   tbox::plog<<"CompositionStrategyMobilities... "<<endl;
   CompositionStrategyMobilities composition_strategy_mobilities(dcalphad_db,
         false,
         1,
         &free_energy_strategy );

   composition_strategy_mobilities.printDiagnostics(temperature,temperature);


   cafe.energyVsPhiAndC(temperature, &lceq[0], found_ceq, with_third_phase, 101, 100);

   input_db.reset();


   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return(0);
}

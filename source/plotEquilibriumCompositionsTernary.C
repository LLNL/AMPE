// Copyright (c) 2018, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and
// the Oak Ridge National Laboratory
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
// LLC, UT BATTELLE, LLC, 
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#include "CALPHADFreeEnergyFunctionsTernary.h"
#include "ConcInterpolationType.h"

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/Database.h"

#include <string>
#include <map>
#include <iostream>
#include <fstream>

using namespace SAMRAI;
using namespace std;


int main( int argc, char *argv[] )
{
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();

   {
   std::string input_filename( argv[1] );

   // Create input database and parse all data in input file.
   boost::shared_ptr<tbox::MemoryDatabase> input_db(
      new tbox::MemoryDatabase("input_db"));
   tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

   // make from input file name
   std::string run_name =
      input_filename.substr( 0, input_filename.rfind( "." ) );

   // Logfile
   std::string log_file_name = run_name + ".log";
   tbox::PIO::logOnlyNodeZero( log_file_name );

#ifdef GITVERSION
#define xstr(x) #x
#define LOG(x) tbox::plog<<" AMPE: git version "<<xstr(x)<<endl;
    LOG(GITVERSION);
    tbox::plog<<endl;
#endif

   tbox::plog << "input_filename = " << input_filename << endl;

   boost::shared_ptr<tbox::Database> model_db =
      input_db->getDatabase("ModelParameters");

   string energy_interp_func_type = "pbg";
   ConcInterpolationType conc_interp_func_type = ConcInterpolationType::PBG;
 
   boost::shared_ptr<tbox::Database> temperature_db =
      model_db->getDatabase( "Temperature" );
   double temperature_low  = temperature_db->getDouble( "low" );
   double temperature_high = temperature_db->getDouble( "high" );

   boost::shared_ptr<tbox::Database> conc_db(
      model_db->getDatabase( "ConcentrationModel" ));

   boost::shared_ptr<tbox::Database> dcalphad_db=
      conc_db->getDatabase( "Calphad" );
   std::string calphad_filename = dcalphad_db->getString( "filename" );
   boost::shared_ptr<tbox::MemoryDatabase> calphad_db (
      new tbox::MemoryDatabase( "calphad_db" ) );
   tbox::InputManager::getManager()->parseInputFile(
      calphad_filename, calphad_db );
   
   boost::shared_ptr<tbox::Database> newton_db;
   if ( conc_db->isDatabase( "NewtonSolver" ) )
      newton_db = conc_db->getDatabase( "NewtonSolver" );

   bool with_third_phase=false;
   
   CALPHADFreeEnergyFunctionsTernary
      cafe(calphad_db, newton_db,
           energy_interp_func_type,
           conc_interp_func_type);

   // choose pair of phases: phaseL, phaseA, phaseB
   const PHASE_INDEX pi0=phaseL;
   const PHASE_INDEX pi1=phaseA;
   
   // initial guesses
   double init_guess[5];
   model_db->getDoubleArray("initial_guess", &init_guess[0], 5);

   double nominalc[2];
   model_db->getDoubleArray("concentration",&nominalc[0],2);
   double lceq[5]={init_guess[0], init_guess[1], // liquid
                   init_guess[2], init_guess[3], // solid
                   init_guess[4]};

   map<double,double> cseq0;
   map<double,double> cseq1;
   map<double,double> cleq0;
   map<double,double> cleq1;

   double dT= (temperature_high-temperature_low)/50;

   //loop over temperature range
   for(int iT=0; iT<50; iT++){

      double temperature = temperature_low + iT*dT;

      // compute equilibrium concentrations
      bool found_ceq =
         cafe.computeCeqT(temperature,pi0,pi1,nominalc[0],
                       nominalc[1],&lceq[0]);
      if( lceq[0]>1. )found_ceq = false;
      if( lceq[0]<0. )found_ceq = false;
      if( lceq[1]>1. )found_ceq = false;
      if( lceq[1]<0. )found_ceq = false;
   
      if( found_ceq ){
         //tbox::pout<<"Found equilibrium concentrations: "
         //          <<lceq[0]<<" and "<<lceq[1]<<"..."<<endl;
         cleq0.insert( pair<double,double>(temperature,lceq[0]) );
         cleq1.insert( pair<double,double>(temperature,lceq[1]) );
         cseq0.insert( pair<double,double>(temperature,lceq[2]) );
         cseq1.insert( pair<double,double>(temperature,lceq[3]) );

      }else{
         tbox::pout<<"Temperature = "<<temperature<<endl;
         tbox::pout<<"ERROR: Equilibrium concentrations not found... "<<endl;
         return 1;
      }

   }

   ofstream os("CvsT.dat");
   os<<"#liquid0\n";
   {
      map<double,double>::iterator it = cleq0.begin();
      while( it!=cleq0.end() ){
         os<<it->first<<"  "<<it->second<<endl;
         ++it;
      }
   }
   os<<endl<<endl;

   os<<"#liquid1\n";
   {
      map<double,double>::iterator it = cleq1.begin();
      while( it!=cleq1.end() ){
         os<<it->first<<"  "<<it->second<<endl;
         ++it;
      }
   }
   os<<endl<<endl;

   os<<"#solid0\n";
   {
      map<double,double>::iterator it = cseq0.begin();
      while( it!=cseq0.end() ){
         os<<it->first<<"  "<<it->second<<endl;
         ++it;
      }
   }
   os<<endl<<endl;

   os<<"#solid1\n";
   {
      map<double,double>::iterator it = cseq1.begin();
      while( it!=cseq1.end() ){
         os<<it->first<<"  "<<it->second<<endl;
         ++it;
      }
   }

   input_db.reset();

   }

   tbox::SAMRAIManager::shutdown();
   tbox::SAMRAIManager::finalize();
   tbox::SAMRAI_MPI::finalize();

   return 0;
}

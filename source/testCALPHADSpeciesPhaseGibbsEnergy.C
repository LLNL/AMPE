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
#include "CALPHADSpeciesPhaseGibbsEnergy.h"

#include "SAMRAI/tbox/MemoryDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/SAMRAIManager.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace SAMRAI;

int main( int argc, char *argv[] )
{
   cout<<"Test CALPHAD Species Gibbs Energy functions."<<endl;

   /*
    * Initialize MPI, SAMRAI, and enable logging.
    */
   tbox::SAMRAI_MPI::init(&argc, &argv);
   tbox::SAMRAIManager::initialize();
   tbox::SAMRAIManager::startup();


   std::string input_filename;
   input_filename = argv[1];
   cout<<"Input file = "<<input_filename<<endl;

   boost::shared_ptr<tbox::MemoryDatabase> input_db =
      tbox::InputManager::getManager()->parseInputFile(input_filename);
   input_db->printClassData(cout);

   if ( input_db->keyExists( "ModelParameters" ) ){
      boost::shared_ptr<tbox::Database> model_db =
         input_db->getDatabase("ModelParameters");

      cout<<"Read T"<<endl;
      boost::shared_ptr<tbox::Database> temperature_db = model_db->getDatabase( "Temperature" );

      double temperature = temperature_db->getDouble( "temperature" );
      cout<<"T="<<temperature<<endl;

      boost::shared_ptr<tbox::Database> conc_db(model_db->getDatabase( "ConcentrationModel" ));
      boost::shared_ptr<tbox::Database> dcalphad_db=conc_db->getDatabase( "Calphad" );
      std::string calphad_filename = dcalphad_db->getString( "filename" );
      boost::shared_ptr<tbox::MemoryDatabase> calphad_db ( new tbox::MemoryDatabase( "calphad_db" ) );
      tbox::InputManager::getManager()->parseInputFile( calphad_filename, calphad_db );

      boost::shared_ptr<tbox::Database> speciesC_db = calphad_db->getDatabase( "SpeciesC" );
      speciesC_db->printClassData(cout);

      string name = speciesC_db->getStringWithDefault( "name", "unknown" );
      
      CALPHADSpeciesPhaseGibbsEnergy gspecies;
      gspecies.initialize(name,speciesC_db->getDatabase( "PhaseL" ) );

      double energy = gspecies.fenergy(temperature);
      cout<<"Energy = "<<energy<<endl;
   }else{
      cout<<"No ModelParameters found in input file"<<endl;
      return 1;
  }

   return(0);
}

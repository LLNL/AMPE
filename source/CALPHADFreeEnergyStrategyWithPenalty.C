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
#include "CALPHADFreeEnergyStrategyWithPenalty.h";
#include "CALPHADConcSolverWithPenalty.h"
#include "CALPHADFreeEnergyFunctionsWithPenaltyBinary.h"

#include "SAMRAI/tbox/InputManager.h"

#include <string>
using namespace std;

using namespace SAMRAI;

CALPHADFreeEnergyStrategyWithPenalty::CALPHADFreeEnergyStrategyWithPenalty(
   boost::shared_ptr<tbox::Database> calphad_db,
   boost::shared_ptr<tbox::Database> newton_db,
   const std::string& phase_interp_func_type,
   const std::string& eta_interp_func_type,
   const std::string& avg_func_type,
   MolarVolumeStrategy* mvstrategy,
   const int conc_l_id,
   const int conc_a_id,
   const int conc_b_id,
   const bool with_third_phase,
   const double  phase_well_scale,
   const double eta_well_scale,
   const std::string& phase_well_func_type,
   const std::string& eta_well_func_type ):
      CALPHADFreeEnergyStrategy(calphad_db,newton_db,
         phase_interp_func_type,
         eta_interp_func_type,
         avg_func_type,
         mvstrategy,
         conc_l_id,
         conc_a_id,
         conc_b_id,
         with_third_phase,
         phase_well_scale,
         eta_well_scale,
         phase_well_func_type,
         eta_well_func_type)
{
   tbox::pout << "CALPHADFreeEnergyStrategyWithPenalty()..." << endl;
   const short n = with_third_phase ? 3 : 2;
   d_penalty_parameters.resize( n );
   for(short i=0;i<n;i++)d_penalty_parameters[i].resize(6);
   
   std::string namemixL("PenaltyPhaseL");
   boost::shared_ptr<tbox::Database> mixL_db = calphad_db->getDatabase( namemixL );
   mixL_db->getDoubleArray( "Left",  &d_penalty_parameters[0][0], 3 );
   mixL_db->getDoubleArray( "Right", &d_penalty_parameters[0][3], 3 );

   std::string namemixA("PenaltyPhaseA");
   boost::shared_ptr<tbox::Database> mixA_db = calphad_db->getDatabase( namemixA );
   mixA_db->getDoubleArray( "Left",  &d_penalty_parameters[1][0], 3 );
   mixA_db->getDoubleArray( "Right", &d_penalty_parameters[1][3], 3 );

   if ( with_third_phase ) {
      std::string namemixB("PenaltyPhaseB");
      boost::shared_ptr<tbox::Database> mixB_db = calphad_db->getDatabase( namemixB );
      mixB_db->getDoubleArray( "Left",  &d_penalty_parameters[2][0], 3 );
      mixB_db->getDoubleArray( "Right", &d_penalty_parameters[2][3], 3 );
   }
   
   setup(calphad_db,newton_db);
}

//=======================================================================

void CALPHADFreeEnergyStrategyWithPenalty::setup(boost::shared_ptr<tbox::Database> calphad_db,
                                                 boost::shared_ptr<tbox::Database> newton_db)
{
   tbox::pout << "CALPHADFreeEnergyStrategyWithPenalty::setupSolver()..." << endl;
   if(d_calphad_fenergy!=0)delete d_calphad_fenergy;
   
   d_calphad_fenergy = new
      CALPHADFreeEnergyFunctionsWithPenaltyBinary(calphad_db,newton_db,d_phase_interp_func_type,
                                 d_eta_interp_func_type,d_avg_func_type,
                                 d_with_third_phase,
                                 d_phase_well_scale,d_eta_well_scale,
                                 d_phase_well_func_type,d_eta_well_func_type);
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature
bool CALPHADFreeEnergyStrategyWithPenalty::computeCeqT(
   const double temperature,
   const PHASE_INDEX pi0, const PHASE_INDEX pi1,
   double* ceq )
{
   assert( temperature>0. );

   return d_calphad_fenergy->computeCeqT(temperature,pi0,pi1,ceq);
}

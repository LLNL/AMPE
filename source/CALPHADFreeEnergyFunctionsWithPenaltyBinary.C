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
#include "CALPHADFreeEnergyFunctionsWithPenaltyBinary.h"
#include "CALPHADConcSolverBinaryWithPenalty.h"
#include "CALPHADEqConcSolverBinaryWithPenalty.h"
#include "PhysicalConstants.h"

using namespace std;

CALPHADFreeEnergyFunctionsWithPenaltyBinary::CALPHADFreeEnergyFunctionsWithPenaltyBinary(
   boost::shared_ptr<SAMRAI::tbox::Database> calphad_db,
   boost::shared_ptr<SAMRAI::tbox::Database> newton_db,
   const std::string& phase_interp_func_type,
   const std::string& eta_interp_func_type,
   const std::string& avg_func_type,
   const bool with_third_phase,
   const double  phase_well_scale,
   const double eta_well_scale,
   const std::string& phase_well_func_type,
   const std::string& eta_well_func_type ):
      CALPHADFreeEnergyFunctionsBinary(
         calphad_db,newton_db,
         phase_interp_func_type,
         eta_interp_func_type,
         avg_func_type,
         with_third_phase,
         phase_well_scale,
         eta_well_scale,
         phase_well_func_type,
         eta_well_func_type)
{
   const short n = with_third_phase ? 3 : 2;
   d_penalty_parameters.resize( n );
   for(short i=0;i<n;i++)d_penalty_parameters[i].resize(6);
   
   readParameters(calphad_db);

   setupSolver(newton_db);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsWithPenaltyBinary::setupSolver(boost::shared_ptr<tbox::Database> newton_db)
{
   tbox::pout << "CALPHADFreeEnergyFunctionsWithPenaltyBinary::setupSolver()..." << std::endl;
   
   if( d_solver!=NULL )delete d_solver;
   d_solver = new CALPHADConcentrationSolverBinaryWithPenalty( d_with_third_phase,
                                                         d_penalty_parameters );

   readNewtonparameters(newton_db);
}

//=======================================================================

void CALPHADFreeEnergyFunctionsWithPenaltyBinary::readParameters(
   boost::shared_ptr<SAMRAI::tbox::Database> calphad_db)
{
   std::string namemixL("PenaltyPhaseL");
   boost::shared_ptr<tbox::Database> mixL_db = calphad_db->getDatabase( namemixL );
   mixL_db->getDoubleArray( "Left",  &d_penalty_parameters[0][0], 3 );
   mixL_db->getDoubleArray( "Right", &d_penalty_parameters[0][3], 3 );

   std::string namemixA("PenaltyPhaseA");
   boost::shared_ptr<tbox::Database> mixA_db = calphad_db->getDatabase( namemixA );
   mixA_db->getDoubleArray( "Left",  &d_penalty_parameters[1][0], 3 );
   mixA_db->getDoubleArray( "Right", &d_penalty_parameters[1][3], 3 );

   if ( d_with_third_phase ) {
      std::string namemixB("PenaltyPhaseB");
      boost::shared_ptr<tbox::Database> mixB_db = calphad_db->getDatabase( namemixB );
      mixB_db->getDoubleArray( "Left",  &d_penalty_parameters[2][0], 3 );
      mixB_db->getDoubleArray( "Right", &d_penalty_parameters[2][3], 3 );
   }

}

//-----------------------------------------------------------------------

double CALPHADFreeEnergyFunctionsWithPenaltyBinary::computeFreeEnergy(
   const double temperature,
   const double* const conc,
   const PHASE_INDEX pi,
   const bool gp )
{
   double fe = CALPHADFreeEnergyFunctionsBinary::computeFreeEnergy(temperature,
      conc,pi,false);

   double extra_energy = computePenalty(pi,conc[0]);
   
   // subtract -mu*c to get grand potential
   if( gp )fe -= computeDerivFreeEnergy(temperature,conc[0],pi)*conc[0];
   
   return fe+extra_energy;
}

//=======================================================================

double CALPHADFreeEnergyFunctionsWithPenaltyBinary::computeDerivFreeEnergy(
   const double temperature,
   const double conc,
   const PHASE_INDEX pi )
{
   double fe = CALPHADFreeEnergyFunctionsBinary::computeDerivFreeEnergy(temperature,
      conc,pi);

   double extra_energy = computeDerivPenalty(pi,conc);
   
   return fe+extra_energy;
}

//=======================================================================

void CALPHADFreeEnergyFunctionsWithPenaltyBinary::computeSecondDerivativeFreeEnergy(
   const double temp,
   const std::vector<double>& conc,
   const PHASE_INDEX pi,
   std::vector<double>& d2fdc2)
{
   CALPHADFreeEnergyFunctionsBinary::computeSecondDerivativeFreeEnergy(temp,
      conc,pi,d2fdc2);

   double extra_energy = compute2ndDerivPenalty(pi, conc[0]);
      
   d2fdc2[0]+= extra_energy;
}

//=======================================================================

// compute equilibrium concentrations in various phases for given temperature
bool CALPHADFreeEnergyFunctionsWithPenaltyBinary::computeCeqT(
   const double temperature,
   const PHASE_INDEX pi0, const PHASE_INDEX pi1,
   double* ceq )
{
   assert( temperature>0. );
   assert( d_penalty_parameters.size()>0 );

   int N = 2;

   double* fA = new double[N];
   double* fB = new double[N];
   double* L0 = new double[N];
   double* L1 = new double[N];
   double* L2 = new double[N];
   double* L3 = new double[N];

   setupValuesForTwoPhasesSolver(temperature, L0, L1, L2, L3, fA, fB, pi0, pi1);
   vector<vector<double> > penalty_parameters;
   penalty_parameters.push_back( d_penalty_parameters[pi0] );
   penalty_parameters.push_back( d_penalty_parameters[pi1] );
   double RTinv = 1.0 / ( gas_constant_R_JpKpmol * temperature );
   CALPHADEqConcentrationSolverBinaryWithPenalty  eq_solver;
   int ret = eq_solver.ComputeConcentrationWithPenalty(
      ceq,
      RTinv,
      L0, L1, L2, L3,
      fA, fB,
      penalty_parameters );

   if( ret>=0 )
   {
      switch( pi0 )
      {
         case phaseL:
            tbox::pout<<"CALPHAD with Penalty, c_eq phaseL="<<ceq[0]<<endl;
            d_ceq_l=ceq[0];
            break;
         case phaseA:
            tbox::pout<<"CALPHAD with Penalty, c_eq phaseA="<<ceq[0]<<endl;
            d_ceq_a=ceq[0];
            break;
         case phaseB:
            tbox::pout<<"CALPHAD with Penalty, c_eq phaseB="<<ceq[0]<<endl;
            d_ceq_b=ceq[0];
            break;
      }
      
      switch( pi1 )
      {
         case phaseL:
            tbox::pout<<"CALPHAD with Penalty, c_eq phaseL="<<ceq[1]<<endl;
            d_ceq_l=ceq[1];
            break;
         case phaseA:
            tbox::pout<<"CALPHAD with Penalty, c_eq phaseA="<<ceq[1]<<endl;
            d_ceq_a=ceq[1];
            break;
         case phaseB:
            tbox::pout<<"CALPHAD with Penalty, c_eq phaseB="<<ceq[1]<<endl;
            d_ceq_b=ceq[1];
            break;
      }
   }
   
   delete[] fA;
   delete[] fB;
   delete[] L0;
   delete[] L1;
   delete[] L2;
   delete[] L3;
   
   return (ret>=0);
}

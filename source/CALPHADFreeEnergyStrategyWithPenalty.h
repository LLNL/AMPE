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
#ifndef included_CALPHADFreeEnergyStrategyWithPenalty
#define included_CALPHADFreeEnergyStrategyWithPenalty

#include "CALPHADFreeEnergyStrategyBinary.h"
#include "CALPHADFunctions.h"
#include "MolarVolumeStrategy.h"


class CALPHADFreeEnergyStrategyWithPenalty:
   public CALPHADFreeEnergyStrategyBinary
{
public:
   CALPHADFreeEnergyStrategyWithPenalty(
      boost::shared_ptr<tbox::Database> input_db,
      boost::shared_ptr<tbox::Database> newton_db,
      const std::string& phase_interp_func_type,
      const std::string& eta_interp_func_type,
      const std::string& avg_func_type,
      MolarVolumeStrategy* mvstrategy,
      const int conc_l_id,
      const int conc_a_id,
      const int conc_b_id,
      const int ncompositions,
      const bool with_third_phase,
      const double  phase_well_scale,
      const double eta_well_scale,
      const std::string& phase_well_func_type,
      const std::string& eta_well_func_type );
      
   ~CALPHADFreeEnergyStrategyWithPenalty(){};
 
   virtual void setup(boost::shared_ptr<tbox::Database> input_db,
                      boost::shared_ptr<tbox::Database> newton_db);

   virtual bool computeCeqT(
      const double temperature,
      const PHASE_INDEX pi0, const PHASE_INDEX pi1,
      double* ceq );

   virtual double computeValFreeEnergyLiquid(
      const double temperature,
      const double conc,
      const bool gp = false )
   {
      const double f1 = d_calphad_fenergy->computeFreeEnergy(temperature,&conc,phaseL,gp);
      const double f2 = d_calphad_fenergy->computePenalty(phaseL, conc);
      
      return (f1+f2)*d_mv_strategy->computeInvMolarVolume(temperature,&conc,phaseL); 
   }

   virtual double computeValFreeEnergySolidA(
      const double temperature,
      const double conc,
      const bool gp = false )
   {
      const double f1 = d_calphad_fenergy->computeFreeEnergy(temperature,&conc,phaseA,gp);
      const double f2 = d_calphad_fenergy->computePenalty(phaseA, conc);
      
      return (f1+f2)*d_mv_strategy->computeInvMolarVolume(temperature,&conc,phaseA); 
   }

   virtual double computeValFreeEnergySolidB(
      const double temperature,
      const double conc,
      const bool gp = false )
   {
      const double f1 = d_calphad_fenergy->computeFreeEnergy(temperature,&conc,phaseB,gp);
      const double f2 = d_calphad_fenergy->computePenalty(phaseB, conc);
      
      return (f1+f2)*d_mv_strategy->computeInvMolarVolume(temperature,&conc,phaseB);  
   }

   void computeSecondDerivativeEnergyPhaseL(
      const double temperature,
      const std::vector<double>& c_l,
      std::vector<double>& d2fdc2,
      const bool use_internal_units)
   {
      CALPHADFreeEnergyFunctionsBinary* calphad_fenergy=dynamic_cast<CALPHADFreeEnergyFunctionsBinary*>(d_calphad_fenergy);
      assert( calphad_fenergy );
      
      calphad_fenergy->computeSecondDerivativeFreeEnergy(temperature,&c_l[0],phaseL,d2fdc2);
      
      double extra_energy = calphad_fenergy->compute2ndDerivPenalty(phaseL, c_l[0]);
      
      d2fdc2[0]+= extra_energy;
  
      if( use_internal_units )
         d2fdc2[0] *= d_mv_strategy->computeInvMolarVolume(temperature,&c_l[0],phaseL);
      
      if( d2fdc2[0]<0. )
         tbox::pout<<"CALPHADFreeEnergyStrategyWithPenalty --- WARNING: fcc="<<d2fdc2[0]<<" in liquid for cl="<<c_l[0]<<"!!!"<<std::endl;
      //tbox::pout<<"fcc="<<d2fdc2[0]<<" in liquid!!!"<<endl;
   }

   void computeSecondDerivativeEnergyPhaseA(
      const double temperature,
      const std::vector<double>& c_a,
      std::vector<double>& d2fdc2,
      const bool use_internal_units)
   {
      CALPHADFreeEnergyFunctionsBinary* calphad_fenergy=dynamic_cast<CALPHADFreeEnergyFunctionsBinary*>(d_calphad_fenergy);
      
      calphad_fenergy->computeSecondDerivativeFreeEnergy(temperature,&c_a[0],phaseA,d2fdc2);
      
      double extra_energy = calphad_fenergy->compute2ndDerivPenalty(phaseA, c_a[0]);
      
      d2fdc2[0]+= extra_energy;
      
      if( use_internal_units )
         d2fdc2[0] *= d_mv_strategy->computeInvMolarVolume(temperature,&c_a[0],phaseA);
      
      if( d2fdc2[0]<0. )
         tbox::pout<<"CALPHADFreeEnergyStrategyWithPenalty --- WARNING: fcc="<<d2fdc2[0]<<" in phase A for ca="<<c_a[0]<<"!!!"<<std::endl;
      //tbox::pout<<"fcc="<<d2fdc2[0]<<" in liquid!!!"<<endl;
   }

   void computeSecondDerivativeEnergyPhaseB(
      const double temperature,
      const std::vector<double>& c_b,
      std::vector<double>& d2fdc2,
      const bool use_internal_units)
   {
      CALPHADFreeEnergyFunctionsBinary* calphad_fenergy=dynamic_cast<CALPHADFreeEnergyFunctionsBinary*>(d_calphad_fenergy);
      
      calphad_fenergy->computeSecondDerivativeFreeEnergy(temperature,&c_b[0],phaseB,d2fdc2);
      
      double extra_energy = calphad_fenergy->compute2ndDerivPenalty(phaseB, c_b[0]);
      
      d2fdc2[0]+= extra_energy;
      
      if( use_internal_units )
         d2fdc2[0] *= d_mv_strategy->computeInvMolarVolume(temperature,&c_b[0],phaseB);
      
      if( d2fdc2[0]<0. )
         tbox::pout<<"CALPHADFreeEnergyStrategyWithPenalty --- WARNING: fcc="<<d2fdc2[0]<<" in phase B for cb="<<c_b[0]<<"!!!"<<std::endl;
      //tbox::pout<<"fcc="<<d2fdc2[0]<<" in liquid!!!"<<endl;
   }

private:

   std::vector< std::vector<double> > d_penalty_parameters;
   
};

#endif


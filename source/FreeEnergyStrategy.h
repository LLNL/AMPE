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
#ifndef included_FreeEnergyStrategy
#define included_FreeEnergyStrategy

#include "PhysicalConstants.h"
#include "Phases.h"

#include "SAMRAI/hier/PatchHierarchy.h"

#include <boost/make_shared.hpp>
#include <vector>

using namespace SAMRAI;


class FreeEnergyStrategy
{
public:
   FreeEnergyStrategy(){};
   virtual ~FreeEnergyStrategy(){};

   // mesh functions
   virtual void computeFreeEnergyLiquid(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int f_l_id,
      const bool gp ) = 0;

   virtual void computeDerivFreeEnergyLiquid(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int f_l_id )
   {
      (void) temperature_id;  // unused
      computeDerivFreeEnergy(hierarchy,f_l_id);
   }

   virtual void computeFreeEnergySolidA(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int f_a_id,
      const bool gp ) = 0;

   virtual void computeDerivFreeEnergySolidA(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int f_a_id )
   {
      (void) temperature_id;  // unused
      computeDerivFreeEnergy(hierarchy,f_a_id);
   }

   virtual void computeFreeEnergySolidB(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int f_b_id,
      const bool gp ) = 0;

   // mesh functions
   virtual void computeFreeEnergyLiquid(
      hier::Patch& patch,
      const int temperature_id,
      const int f_l_id,
      const bool gp ) = 0;

   virtual void computeFreeEnergySolidA(
      hier::Patch& patch,
      const int temperature_id,
      const int f_a_id,
      const bool gp ) = 0;

   virtual void computeFreeEnergySolidB(
      hier::Patch& patch,
      const int temperature_id,
      const int f_b_id,
      const bool gp ) = 0;

   virtual void addComponentRhsPhi(
      hier::Patch& patch,
      const int temperature_id,
      const int phase_id,
      const int eta_id,
      const int conc_id, 
      const int f_l_id,
      const int f_a_id,
      const int f_b_id,
      const int rhs_id ) = 0;

   virtual void addComponentRhsEta(
      hier::Patch& patch,
      const int temperature_id,
      const int phase_id,
      const int eta_id,
      const int conc_id, 
      const int f_l_id,
      const int f_a_id,
      const int f_b_id,
      const int rhs_id ) = 0;

   // pointwise functions
   virtual void computeSecondDerivativeEnergyPhaseL(
      const double temperature,
      const std::vector<double>& c,
      std::vector<double>& d2fdc2, const bool use_internal_units=true) = 0;
   virtual void computeSecondDerivativeEnergyPhaseA(
      const double temperature,
      const std::vector<double>& c,
      std::vector<double>& d2fdc2, const bool use_internal_units=true) = 0;
   virtual void computeSecondDerivativeEnergyPhaseB(
      const double temperature,
      const std::vector<double>& c,
      std::vector<double>& d2fdc2, const bool use_internal_units=true) = 0;

   virtual void preRunDiagnostics(std::ostream& os)
   {
      return;
   }

private:

   void computeDerivFreeEnergy(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int df_id );


};
#endif

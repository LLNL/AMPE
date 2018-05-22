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
#ifndef included_BiasDoubleWellFreeEnergyStrategy
#define included_BiasDoubleWellFreeEnergyStrategy 

#include "FreeEnergyStrategy.h"
#include "MeltingTemperatureStrategy.h"

/*!
 * @brief Class BiasDoubleWellFreeEnergyStrategy implements the free energy
 * used by Kobayashi in Physica D 63 (1993), p. 410-423
 * with
 *    m(T)=(alpha/pi)tan^-1(gamma(Te-T))
 *
 */

class BiasDoubleWellFreeEnergyStrategy:
   public FreeEnergyStrategy
{
public:
   BiasDoubleWellFreeEnergyStrategy(
      const double alpha,
      const double gamma,
      MeltingTemperatureStrategy* meltingTstrat );

   ~BiasDoubleWellFreeEnergyStrategy(){};
 
   void computeFreeEnergyLiquid(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int fl_id,
      const bool gp ){};

   void computeFreeEnergySolidA(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int fs_id,
      const bool gp ){};

   void computeFreeEnergySolidB(
      const boost::shared_ptr<hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int fs_id,
      const bool gp ){};

   // mesh functions
   virtual void computeFreeEnergyLiquid(
      hier::Patch& patch,
      const int temperature_id,
      const int f_l_id,
      const bool gp ){};

   virtual void computeFreeEnergySolidA(
      hier::Patch& patch,
      const int temperature_id,
      const int f_a_id,
      const bool gp ){};

   virtual void computeFreeEnergySolidB(
      hier::Patch& patch,
      const int temperature_id,
      const int f_b_id,
      const bool gp ){};

   void addComponentRhsPhi(
      hier::Patch& patch,
      const int temperature_id,
      const int phase_id,
      const int eta_id,
      const int conc_id, 
      const int fl_id,
      const int fa_id,
      const int fb_id,
      const int rhs_id );

   void addComponentRhsEta(
      hier::Patch& patch,
      const int temperature_id,
      const int phase_id,
      const int eta_id,
      const int conc_id, 
      const int fl_id,
      const int fa_id,
      const int fb_id,
      const int rhs_id );

   void computePhaseConcentrations(
      const boost::shared_ptr< hier::PatchHierarchy > hierarchy,
      const int temperature_id,
      const int phase_id,
      const int eta_id,
      const int concentration_id );

   void computeSecondDerivativeEnergyPhaseL(
      const double temperature,
      const std::vector<double>& c,
      std::vector<double>& d2fdc2, const bool use_internal_units=true);
   void computeSecondDerivativeEnergyPhaseA(
      const double temperature,
      const std::vector<double>& c,
      std::vector<double>& d2fdc2, const bool use_internal_units=true);
   void computeSecondDerivativeEnergyPhaseB(
      const double temperature,
      const std::vector<double>& c,
      std::vector<double>& d2fdc2, const bool use_internal_units=true);

private:

   double d_scale;
   double d_alpha;
   double d_gamma;
   MeltingTemperatureStrategy* d_meltingTstrat;
};

#endif

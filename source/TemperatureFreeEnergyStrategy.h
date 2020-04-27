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
#ifndef included_TemperatureFreeEnergyStrategy
#define included_TemperatureFreeEnergyStrategy

#include "FreeEnergyStrategy.h"
#include "InterpolationType.h"

#include <string>
#include <vector>

class TemperatureFreeEnergyStrategy : public FreeEnergyStrategy
{
 public:
   TemperatureFreeEnergyStrategy(
       const EnergyInterpolationType phase_interp_func_type,
       const EnergyInterpolationType eta_interp_func_type, const double fa,
       const double fb, const double vma, const double vmb,
       const double latent_heat, const double meltingT,
       const bool with_third_phase);

   ~TemperatureFreeEnergyStrategy(){};

   void computeFreeEnergyLiquid(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int fl_id, const bool gp = false);

   void computeFreeEnergySolidA(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int fs_id, const bool gp = false);

   void computeFreeEnergySolidB(
       const std::shared_ptr<hier::PatchHierarchy> hierarchy,
       const int temperature_id, const int fs_id, const bool gp = false);

   void computeFreeEnergyLiquid(hier::Patch& patch, const int temperature_id,
                                const int fl_id, const bool gp = false);

   void computeFreeEnergySolidA(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp = false);

   void computeFreeEnergySolidB(hier::Patch& patch, const int temperature_id,
                                const int fs_id, const bool gp = false);

   void addDrivingForce(const double time, hier::Patch& patch,
                        const int temperature_id, const int phase_id,
                        const int eta_id, const int conc_id, const int fl_id,
                        const int fa_id, const int fb_id, const int rhs_id);

   void addDrivingForceEta(const double time, hier::Patch& patch,
                           const int temperature_id, const int phase_id,
                           const int eta_id, const int conc_id, const int fl_id,
                           const int fa_id, const int fb_id, const int rhs_id);

   void computeSecondDerivativeEnergyPhaseL(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   void computeSecondDerivativeEnergyPhaseA(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);
   void computeSecondDerivativeEnergyPhaseB(
       const double temperature, const std::vector<double>& c,
       std::vector<double>& d2fdc2, const bool use_internal_units = true);

 private:
   double d_latent_heat;
   double d_meltingT;
   double d_invMeltingT;
   double d_f_a;
   double d_f_b;

   EnergyInterpolationType d_phase_interp_func_type;
   EnergyInterpolationType d_eta_interp_func_type;

   bool d_with_third_phase;

   double computeValFreeEnergyLiquid(const double temperature,
                                     const double conc, const bool gp = false);
};

#endif

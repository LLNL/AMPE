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
#ifndef included_PhaseFACOps
#define included_PhaseFACOps

#include "EllipticFACOps.h"
#include "InterpolationType.h"

#include <string>

using namespace SAMRAI;

class PhaseFACOps : public EllipticFACOps
{

 public:
   PhaseFACOps(const std::string& object_name, const bool with_third_phase,
               const std::shared_ptr<tbox::Database> database =
                   std::shared_ptr<tbox::Database>());

   void setOperatorCoefficients(
       const int phase_id, const int eta_id, const int phase_mobility_id,
       const double epsilon_phase, const double gamma,
       const EnergyInterpolationType phase_interp_func_type,
       const double phase_well_scale, const std::string phase_well_func_type,
       const double eta_well_scale = 0.,
       const std::string eta_well_func_type = "");

   void setC(const int phi_id, const int eta_id, const double gamma,
             const EnergyInterpolationType phi_interp_func_type,
             const double phi_well_scale, const std::string phi_well_func_type,
             const double eta_well_scale, const std::string eta_well_func_type);

 private:
   // Timers
   std::shared_ptr<tbox::Timer> t_setcoeffs_timer;

   void setCOnPatchPrivate(std::shared_ptr<pdat::CellData<double> > cd_phi,
                           std::shared_ptr<pdat::CellData<double> > cd_eta,
                           std::shared_ptr<pdat::CellData<double> > cd_m,
                           std::shared_ptr<pdat::CellData<double> > cd_c,
                           const double gamma,
                           const EnergyInterpolationType phi_interp_func_type,
                           const double phi_well_scale,
                           const char* phi_well_func_type,
                           const double eta_well_scale,
                           const char* eta_well_func_type,
                           const hier::Box& pbox);

   bool d_with_third_phase;
};

#endif  // included_PhaseFACOps
